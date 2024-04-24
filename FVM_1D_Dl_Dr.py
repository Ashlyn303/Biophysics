# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 20:52:06 2022


FVM_1D_Sim_Guest_Diffusion_into_Nanopores

@author: Ashlyn

"""
import csv
import math
import numpy as np
from scipy import integrate
from FVM_Function import FVM_Matrices

class FVM_1D_Sim_Dslither_Dl_Dr:

    def __init__(self, Ns, Nz, km, K, ka, kd, Dpore_l, Dpore_r, Dslither, mu, rho, Q_r, Qt, Rg, Lp, Rp0,
                 GAM_0, C0, v0, Bmax, t_array, t_eval, Boundary_Condition, Dirichlet_Conc,
                 Neumann_Flux, Robin_km, Robin_Conc, Labels_X, Labels_Y, Title, color):
        'FVM elements:'
        self.Nz = Nz            # Number of elements along the z-direction

        'Species:'
        self.Ns = Ns            # Number of sepcies is diffusing in the nanopores

        'Fluid:'
        # Mass transfer coefficient, part of the robin boundary condition  (cm/s)
        self.km = km
        # Langmuir adsorption equilibrium constant K = ka/kd   (cm**3/mol)
        self.K = K
        self.ka = ka               # Adsorption coefficient (cm**3/mol s)
        self.kd = kd               # Desorption coefficient (s**-1)
        # Diffusion coeff. of left  side of crystal (cm**2/s)
        self.Dpore_l = Dpore_l
        # Diffusion coeff. of right side of crystal (cm**2/s)
        self.Dpore_r = Dpore_r
        self.Dslither = Dslither   # Diffusion coeff. (cm**2/s)
        self.mu = mu               # Diffusion (viscosity) constant (g/cm/s)
        self.rho = rho             # Density (g/cm**3)
        self.Q_r = Q_r             # Volumetric flow (cm**3/s)
        # Time varying part of volumetric flow (no units)
        self.Qt = Qt

        'Geometry:'
        self.Rg = Rg             # Guest radius (cm)
        self.Lp = Lp             # Length in x-direction (cm)
        self.Rp0 = Rp0           # Length in r-direction (cm)

        'Forces:'
        ''' Array Input: GAM_0 '''
        self.GAM_0 = GAM_0       # Velocity-independent force per volume (g*cm/s**2)/cm**3

        'Initial and Boundary Conditions:'
        ''' Array Input: C0, v0 '''
        self.C0 = C0             # Initial concentration
        self.v0 = v0             # Initial velocity in z-direction (cm/s)
        self.Bmax = Bmax         # Maximum adsorbed conc.(mol/cm^3)

        'Evalution time:'
        self.t_array = t_array   # Time span
        self.t_eval = t_eval     # Time points

        'Type of Boundary Condition:'
        # Set the type of boundary condition
        self.Boundary_Condition = Boundary_Condition
        'Dirichlet Boundary Condition:'
        self.Dirichlet_Conc = Dirichlet_Conc     # Dirichlet_Conc = Csol
        'Neumann Boundary Condition:'
        self.Neumann_Flux = Neumann_Flux         # Value of fixed fluxes J*A = mg/s
        'Robin Boundary Condition:'
        # Mass transfer rate for these Robin BC's (km*A = cm**3/s)
        self.Robin_km = Robin_km
        # External concentrations at these Robin BC's (mg/cm**3)
        self.Robin_Conc = Robin_Conc

        'Plotting properties:'
        'Create labels and title for plot'
        self.Labels_X = Labels_X
        self.Labels_Y = Labels_Y
        self.Title = Title
        self.color = color

    def parameters(self):
        'Get discritization points (edges)'
        x = np.linspace(0, self.Lp, self.Nz+1, endpoint=True)
        self.xarray = np.tile(x, self.Ns)
        'Get element center points'
        'Compute the location of the midpoints for each element.'
        xmiddle = (self.xarray[1:self.Nz + 1] + self.xarray[0:self.Nz])/2
        self.xmid = np.tile(xmiddle, self.Ns)

        'Assign linear adsorption reaction kinetic coeffiction ka to FVM transport reaction term GAM_1A (A -(ka)-> B, bounding)'
        self.GAM_1A = self.ka

        'Assign linear desorption reaction kinetic coeffiction kd to FVM transport reaction term GAM_1B (B --> A, unbounding)'
        self.GAM_1B = self.kd

        'Get initial concentration for each elements at center points'
        # self.C0_mid = self.C0(self.xmid)
        self.C0_mid = self.C0
        return self.xarray, self.xmid, self.GAM_1A, self.GAM_1B, self.C0_mid

    def Global_FVM_Matrices(self, t, C):
        self.parameters()

        'Form FVM Matrices'
        'Allocate FVM matrices & FVM BCs:'
        Accumulation = np.zeros(
            (self.Nz*self.Ns, self.Nz*self.Ns), dtype=float)
        Advection = np.zeros((self.Nz*self.Ns, self.Nz*self.Ns), dtype=float)
        Diffusion = np.zeros((self.Nz*self.Ns, self.Nz*self.Ns), dtype=float)
        adsorption = np.zeros((self.Nz), dtype=float)
        React0 = np.zeros((self.Nz*self.Ns), dtype=float)
        Reaction = np.zeros((self.Nz*self.Ns, self.Nz*self.Ns), dtype=float)
        BC_Diff = np.zeros((self.Nz*self.Ns, 2), dtype=float)
        BC_Diff_q = np.zeros((self.Nz*self.Ns, 2), dtype=float)
        BC_q = np.zeros((self.Nz*self.Ns), dtype=float)
        BC_C = np.zeros((self.Nz*self.Ns), dtype=float)

        'Loop through elements and assemble matrices:'
        'Heights of every element.'
        harray = self.xarray[1:self.Nz+1]-self.xarray[0:self.Nz]
        'Fix for element 101 to prevent it get the wrong hi, which is x_array(Nz-1)'
        'And the correct hi for element 101 should be x_array(Nz)'
        new_x_array = np.append(
            self.xarray[0:self.Nz], self.xarray[self.Nz+2:])
        h = np.tile(harray, self.Ns)

        for i in range(0, self.Nz*self.Ns):
            'Find thickness of current and neighboring elements'
            hi = h[i]
            if i == 0 or i == self.Nz:
                him1 = 0
                hip1 = h[i + 1]
            elif i == self.Nz - 1 or i == self.Nz*self.Ns - 1:
                hip1 = 0
                him1 = h[i - 1]
            else:
                him1 = h[i - 1]
                hip1 = h[i + 1]

            if i in range(0, int(round(self.Nz/2))):
                Ci = C[i]           # Ci = 'CA'-unbounded
                qi = C[i+self.Nz]   # qi = 'CB'-bounded

                Dp_i = self.Dpore_l
                Dp_im1 = self.Dpore_l
                Dp_ip1 = self.Dpore_l

                adsorption[i] = Ci * (self.Bmax-qi)

            elif i in range(int(round(self.Nz/2)), self.Nz):
                Ci = C[i]           # Ci = 'CA'-unbounded
                qi = C[i+self.Nz]   # qi = 'CB'-bounded

                Dp_i = self.Dpore_r
                Dp_im1 = self.Dpore_r
                Dp_ip1 = self.Dpore_r

                adsorption[i] = Ci * (self.Bmax-qi)

            elif i in range(self.Nz, self.Nz*self.Ns):

                Dp_i = self.Dslither
                Dp_im1 = self.Dslither
                Dp_ip1 = self.Dslither

            'Derived FVM_Matrices function from MATLAB to get transport matrices for each element'
            [Accum_Jac, Advec_Jac, Difsn_Jac, GAM1A_Jac, GAM1B_Jac, Gam_0_B] = FVM_Matrices(
                Dp_i, Dp_im1, Dp_ip1, 1, self.GAM_1A, self.GAM_1B, self.Q_r, self.Rp0, hi, him1, hip1, new_x_array[i])

            'Assign element transport matrices in the appropriate locations for global FVM matrices:'
            'Accumulation is mass matrix'
            Accumulation[i, i] = Accum_Jac
            'React0 is for Velocity-independent force per volume'
            React0[i] = Gam_0_B

            'Adsorption and Desorption is chemical reaction, bounding and unbounding.'
            if i in range(0, self.Nz):
                Reaction[i, i] = 0  # -GAM1A_Jac
                Reaction[i, i + self.Nz] = GAM1B_Jac
                Reaction[i + self.Nz, i] = 0  # GAM1A_Jac
                Reaction[i + self.Nz, i + self.Nz] = -GAM1B_Jac

            Adsorption = np.concatenate(
                (-GAM1A_Jac * adsorption[:], GAM1A_Jac * adsorption[:]), axis=0).reshape(self.Nz*self.Ns)

            if i == 0:
                Diffusion[i, i:i + 2] = Difsn_Jac[1:3]
                Dirichelet_bc_l = Difsn_Jac[0]
                Advection[i, i:i + 2] = Advec_Jac[1:3]
            elif i == self.Nz - 1:
                Diffusion[i, i - 1:i+1] = Difsn_Jac[0:2]
                Dirichelet_bc_r = Difsn_Jac[-1]
                Advection[i, i - 1:i+1] = Advec_Jac[0:2]
            elif i == self.Nz:
                Diffusion[i, i:i + 2] = Difsn_Jac[1:3]
                Dirichelet_bc_ql = Difsn_Jac[0]
                Advection[i, i:i + 2] = Advec_Jac[1:3]
            elif i == self.Nz*self.Ns - 1:
                Diffusion[i, i - 1:i+1] = Difsn_Jac[0:2]
                Dirichelet_bc_qr = Difsn_Jac[-1]
                Advection[i, i - 1:i+1] = Advec_Jac[0:2]
            else:
                Diffusion[i, i - 1:i + 2] = Difsn_Jac[0:3]
                Advection[i, i - 1:i + 2] = Advec_Jac[0:3]

        def Force(t): return React0*self.GAM_0(t)
        'Apply Boundary Conditions:'
        'Dirichlet BCs:'
        if self.Boundary_Condition == 'Dirichlet':
            
            qL = C[self.Nz]
            qR = C[-1]

            BC_q[self.Nz] = Dirichelet_bc_ql * qL
            BC_q[-1] = Dirichelet_bc_qr * qR
            
            BC_Diff[0, 0] = Dirichelet_bc_l
            BC_Diff[self.Nz-1, 1] = Dirichelet_bc_r

            Boundary_Cs = lambda t: np.dot(BC_Diff, self.Dirichlet_Conc(t))  

        'Neumann BCs:'
        if self.Boundary_Condition == 'Neumann':
            BC_Diff[0, 0] = Dirichelet_bc_l
            BC_Diff[self.Nz-1, 1] = Dirichelet_bc_r
            def Boundary_Cs(t): return np.dot(BC_Diff, self.Neumann_Flux(t))
            'C'
            Diffusion[0, 0] = Diffusion[0, 0]/2
            Diffusion[self.Nz-1, self.Nz-1] = Diffusion[self.Nz-1, self.Nz-1]/2
            'q'
            Diffusion[self.Nz, self.Nz] = Diffusion[self.Nz, self.Nz]/2
            Diffusion[self.Nz*self.Ns-1, self.Nz*self.Ns -
                      1] = Diffusion[self.Nz*self.Ns-1, self.Nz*self.Ns-1]/2

        'Robin BCs:'
        if self.Boundary_Condition == 'Robin':
            
            area = math.pi*self.Rp0**2

            'C'
            BC_Diff[0, 0] = area*self.Robin_km[0]
            BC_Diff[self.Nz-1, 1] = -area*self.Robin_km[1]

            def Boundary_Cs(t): return np.dot(BC_Diff, self.Robin_Conc(t))

            Diffusion[0, 0] = Diffusion[0, 0]/2 - area*self.Robin_km[0]
            Diffusion[self.Nz-1, self.Nz-1] = Diffusion[self.Nz-1, self.Nz-1]/2 + area*self.Robin_km[1]

            'q'
            BC_Diff_q[self.Nz, 0] = area*self.Robin_km[0]
            BC_Diff_q[self.Nz*self.Ns-1, 1] = -area*self.Robin_km[1]

            Diffusion[self.Nz, self.Nz] = Diffusion[self.Nz, self.Nz] / \
                2 - area*self.Robin_km[0]
            Diffusion[self.Nz*self.Ns-1, self.Nz*self.Ns-1] = Diffusion[self.Nz*self.Ns -
                                                                        1, self.Nz*self.Ns-1]/2 + area*self.Robin_km[1]
        """
        Set up the numerical analysis -- the RHS of the ODE:   
            Note the form of this ODE is: 
            dC/dt = inv(M) * [Jac*C(t) + J(t)]                                                          
            where M is the 'mass matrix', Jac is the jacobian, and
            J(t) are defined by boundary conditions.
        """
        return np.dot(np.linalg.inv(Accumulation), (np.dot((Diffusion + Reaction), C) + Adsorption + Boundary_Cs(t) + BC_C + BC_q))

    def get_Results(self):
        self.parameters()
        # print('Solving FVM...')

        'Solve numerical solution of ODE by using integrate.solve_ivp solver'
        method = 'BDF'      # For non-stiff problem
        # method = 'Radau'  # For non-stiff problem
        # method = 'LSODA'  # For non-stiff problem

        # method = 'RK45'   # For stiff problem
        # method = 'RK23'   # For stiff problem
        # method = 'DOP853' # For stiff problem

        print('Integration method =', method)
        FVM_results = integrate.solve_ivp(self.Global_FVM_Matrices, [
                                          self.t_array[0], self.t_array[-1]], self.C0_mid, t_eval=self.t_eval,  method=method, dense_output=True, rtol=2.5e-14, atol=1e-17)

        return FVM_results
