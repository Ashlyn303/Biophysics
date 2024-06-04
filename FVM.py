# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ashlyn
"""
import time
import numpy as np
import csv
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import make_interp_spline, UnivariateSpline
from FVM_1D_Dl_Dr import FVM_1D_Sim_Dslither_Dl_Dr

'---Image settings-----'
size = 25
plt.rc('axes', titlesize=size)     # fontsize of the axes title
plt.rc('axes', labelsize=size)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=size)    # fontsize of the tick labels
plt.rc('legend', fontsize=size)    # legend fontsize
plt.rc('figure', titlesize=5)      # fontsize of the figure title
mpl.rcParams['figure.figsize'] = 4, 5
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

start_time = time.time()

file_path = 'new_Crystal_csv/'  # File path
file_path_2 = 'new_Crystal_smooth_csv/'  # File path

save_file_path = 'savefig/05_FVM/'


'Crystal1'
# xtal = 'C1'                                # Crystal ID
# Lp = 27.2*10**-4                           # Length of nanopore (μm)
# name = file_path+'Crystal1_fold'           # File path to images csv data
# nameBC_l = file_path+'Crystal1_fold_left'  # File path to left boundary images csv data
# nameBC_r = file_path+'Crystal1_fold_right' # File path to right boundary images csv data
# Dpi_l = [3.0439701656946353e-09, 2.0572602922671336e-09, 2.3957736416120542e-09, 2.27393803936258e-09, 1.9084931190703314e-09, 2.708628747564717e-09, 2.8640443749237903e-09, 3.406312669816462e-09]
# Dpi_r = [6.372882544781162e-09, 3.3139655255964863e-09, 2.161954261464783e-09, 2.0675907424289487e-09, 2.3413439567249504e-09, 3.2544289068622564e-09, 3.4431618813761742e-09, 4.301770232736738e-09]
# xtal_fold_to_run = [1, 5]   #[Best fit, Worst fit] SSD=[7.97e-15, 1.46e-14]


'Crystal2'
xtal = 'C2'                              # Crystal ID
Lp = 20.4*10**-4                         # Length of nanopore (μm)
name = file_path+'Crystal2_fold'         # File path to images csv data
nameBC_l = file_path+'Crystal2_fold_left'     # File path to left boundary images csv data
nameBC_r = file_path+'Crystal2_fold_right'    # File path to right boundary images csv data
Dpi_l = [3.7873724528370184e-09, 3.782379304885544e-09, 2.9057318978610903e-09, 2.8095959555695562e-09, 3.316526603608882e-09, 3.130296230995199e-09, 3.9177986733176005e-09, 4.667938966173201e-09]
Dpi_r = [4.743973417099125e-09, 2.880165676066599e-09, 2.9041696560652437e-09, 2.801303183615771e-09, 2.426148850482564e-09, 1.985368144836361e-09, 2.0655052294953742e-09, 2.808836168379461e-09]
xtal_fold_to_run = [4, 1]  #[Best fit, Worst fit] SSD=[3.62e-15, 6.77e-15]

'Crystal3'
# xtal = 'C3'                              # Crystal ID
# Lp = 61.88*10**-4                        # Length of nanopore (μm)
# name = file_path+'Crystal3_fold'         # File path to images csv data
# nameBC_l = file_path+'Crystal3_fold_left'     # File path to left boundary images csv data
# nameBC_r = file_path+'Crystal3_fold_right'    # File path to right boundary images csv data
# Dpi_l = [1.08388336330001e-08, 1.1392962467143132e-08, 1.0411626663017834e-08, 9.599432805857126e-09, 1.0798327995178489e-08, 1.1918750586006434e-08, 1.4580069399939753e-08, 1.7158081786293558e-08]
# Dpi_r = [1.4328297126743272e-08, 1.1470087790857829e-08, 1.2561026433600456e-08, 1.1998487945197064e-08, 1.13528943434993e-08, 1.2337262533992742e-08, 1.2083042156800448e-08, 1.3100433347135533e-08]
# xtal_fold_to_run = [4, 8]   #[Best fit, Worst fit] SSD=[6.43e-15, 1.60e-14]

print('Average Diffusion Coefficient = {0:.2e}±{1:.2e}'.format(np.mean(Dpi_l+Dpi_r), np.std(Dpi_l+Dpi_r)))

' - - NUMERIC VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - '
Nz = 100
Ns = 2
K = 4.18e9
ka = 0
h = 0.68*1*10*1e-4  # cm/pixel = 6.8e-5
km = 1/h
kd = ka/K
mu = 1e-3
rho = 1
Q_r = 0
Rg = 8.1e-8 
Rp0 = 6.50e-7
def Qt(t): return 1
def GAM_0(t): return 0+0*t
def v0(r): return 0+0*r
Csol = 5e-10
Bmax = 3.13e-5

err_vec = []
# for xtal_fold_i in range(8):
    # xtal_fold_i+=1
for xtal_fold_i in xtal_fold_to_run:
    
    Dpore_l = Dpi_l[xtal_fold_i-1]
    Dpore_r = Dpi_r[xtal_fold_i-1]
    Dslither = 0

    t_array = [0, 1820]

    if xtal == 'C1':
        t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
                  1087, 1207, 1327, 1447, 1567, 1687, 1807]
    elif xtal == 'C2':
        t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
                  1087, 1207, 1327, 1447, 1567, 1687, 1807]

    if xtal == 'C3':
        t_eval = [17, 137, 257, 377, 497, 617, 737, 857, 977,
                  1097, 1217, 1337, 1457, 1577, 1697, 1817]

    pt = [1, 4, 8, 11, -1]

    print('Diffusion rate of left side={0:.2e}'.format(Dpore_l))
    print('Diffusion rate of right side={0:.2e}'.format(Dpore_r))

    "Get BC data from image detect"
    Images_filename = name + str(xtal_fold_i) + '_images_data.csv'
    Images_BC_filename_l = nameBC_l + str(xtal_fold_i) + '_images_BC_data.csv'
    Images_BC_filename_r = nameBC_r + str(xtal_fold_i) + '_images_BC_data.csv'
    # Images_BC_filename_l = nameBC_l + str(xtal_fold_i) + '_images_smooth_BC_data.csv'
    # Images_BC_filename_r = nameBC_r + str(xtal_fold_i) + '_images_smooth_BC_data.csv'
    
    Images_data = pd.read_csv(Images_filename, delimiter=',')

    Images_BC1_data_l = pd.read_csv(Images_BC_filename_l, delimiter=',')
    Image_BC_ave_data_l = Images_BC1_data_l.iloc[:, 0]

    Images_BC1_data_r = pd.read_csv(Images_BC_filename_r, delimiter=',')
    Image_BC_ave_data_r = Images_BC1_data_r.iloc[:, -1][::-1]
    
    C0 = Images_data.values[0, :]
    C0 = np.concatenate((C0, np.repeat(0, Nz)), axis=0)

    'BSpline fit BC'
    Image_BC_Conc_l = Image_BC_ave_data_l
    Image_BC_Conc_r = Image_BC_ave_data_r

    if xtal == 'C1':
        'Left'
        BC_l = Image_BC_Conc_l
        BC_min_l = np.min(BC_l)
        BC_max_l = np.max(BC_l)
        t_time = np.append([0, 0, 0, 0, 0], t_eval)
        t_time = np.append(
            t_time, [t_eval[-1], t_eval[-1], t_eval[-1], t_eval[-1], t_eval[-1]])
        BC_conc_l = np.append(
            [BC_min_l, BC_min_l, BC_min_l, BC_min_l, BC_min_l], BC_l)
        BC_conc_l = np.append(
            BC_conc_l, [BC_max_l, BC_max_l, BC_max_l, BC_max_l, BC_max_l])
        BC_spline_func_l = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_l, size=1), k=3)

        'Right'
        BC_r = Image_BC_Conc_r
        BC_min_r = np.min(BC_r)
        BC_max_r = np.max(BC_r)
        BC_conc_r = np.append(
            [BC_min_r, BC_min_r, BC_min_r, BC_min_r, BC_min_r], BC_r)
        BC_conc_r = np.append(
            BC_conc_r, [BC_max_r, BC_max_r, BC_max_r, BC_max_r, BC_max_r])
        BC_spline_func_r = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_r, size=1), k=3)
        

    elif xtal == 'C2':
        'Left'
        BC_l = Image_BC_Conc_l
        BC_min_l = np.min(BC_l)
        BC_max_l = np.max(BC_l)
        t_time = np.append([0, 0, 0, 0, 0], t_eval)
        BC_conc_l = np.append(
            [BC_min_l, BC_min_l, BC_min_l, BC_min_l, BC_min_l], BC_l)
        BC_spline_func_l = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_l, size=5), k=3)

        'Right'
        BC_r = Image_BC_Conc_r
        BC_min_r = np.min(BC_r)
        BC_max_r = np.max(BC_r)
        BC_conc_r = np.append(
            [BC_min_r, BC_min_r, BC_min_r, BC_min_r, BC_min_r], BC_r)
        BC_spline_func_r = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_r, size=5), k=3)

    elif xtal == 'C3':
        'Left'
        BC_l = Image_BC_Conc_l
        BC_min_l = np.min(BC_l)
        BC_max_l = np.max(BC_l)
        t_time = np.append([0, 0, 0, 0, 0], t_eval)
        t_time = np.append(
            t_time, [t_eval[-1], t_eval[-1], t_eval[-1], t_eval[-1], t_eval[-1]])
        BC_conc_l = np.append(
            [BC_min_l, BC_min_l, BC_min_l, BC_min_l, BC_min_l], BC_l)
        BC_conc_l = np.append(
            BC_conc_l, [BC_max_l, BC_max_l, BC_max_l, BC_max_l, BC_max_l])
        BC_spline_func_l = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_l, size=5), k=3)

        'Right'
        BC_r = Image_BC_Conc_r
        BC_min_r = np.min(BC_r)
        BC_max_r = np.max(BC_r)
        BC_conc_r = np.append(
            [BC_min_r, BC_min_r, BC_min_r, BC_min_r, BC_min_r], BC_r)
        BC_conc_r = np.append(
            BC_conc_r, [BC_max_r, BC_max_r, BC_max_r, BC_max_r, BC_max_r])
        BC_spline_func_r = UnivariateSpline(
            t_time, uniform_filter1d(BC_conc_r, size=5), k=3)

    'Spline function of left boundary and plot'
    t_new = np.linspace(t_eval[0], t_eval[-1], 100)
    c_new = BC_spline_func_l(t_new)
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig1 = plt.figure()
            a0 = fig1.add_axes([0, 0, 1, 1])
            a0.plot(t_eval, uniform_filter1d(Image_BC_Conc_l*1e6, size=1), 'co')
            a0.plot(t_new, c_new*1e6, '-', color='tab:orange', linewidth=3)
            plt.xlabel('Time [sec]')
            plt.ylabel('Concentration [mM]')
            a0.tick_params(direction="in", length=4, width=3)
            fig1.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_FVM_BC_spline_fit_left.svg',
                         format='svg', dpi=400,  bbox_inches='tight')
           
    'Spline function of right boundary and plot'
    t_new = np.linspace(t_eval[0], t_eval[-1], 100)
    c_new = BC_spline_func_r(t_new)
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig1 = plt.figure()
            a0 = fig1.add_axes([0, 0, 1, 1])
            a0.plot(t_eval, uniform_filter1d(Image_BC_Conc_r*1e6, size=1), 'co')
            a0.plot(t_new, c_new*1e6, '-', color='tab:orange', linewidth=3)
            plt.xlabel('Time [sec]')
            plt.ylabel('Concentration [mM]')
            a0.tick_params(direction="in", length=4, width=3)
            fig1.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_FVM_BC_spline_fit_right.svg',
                         format='svg', dpi=400,  bbox_inches='tight')

    'BSpline fit Xtal BC'
    Xtal_BC = Images_data.iloc[:, 0]
    Xtal_BC_spline_func = make_interp_spline(
        t_eval, uniform_filter1d(Xtal_BC, size=1), k=1)

    'Boundary Condition'
    Boundary_Condition = 'Dirichlet'
    # Boundary_Condition = 'Neumann'
    # Boundary_Condition = 'Robin'

    def Dirichlet_Conc(t): return np.array(
        [BC_spline_func_l(t), BC_spline_func_r(t)])
    
    def Neumann_Flux(t): return (BC_spline_func_l(t)-Xtal_BC_spline_func(t))/h

    Robin_km = np.array([km, km])

    # def Robin_Conc(t): return np.array(
    #     [BC_spline_func(t), BC_spline_func(t)])
    def Robin_Conc(t): return np.array(
        [Csol, Csol])
    """
    Input:
    -------
    Create labels and title for plot:
        Labels_X: x-label of figure
        Labels_Y: y-label of figure
        Title: title of figure
        color: color list for plot (Note: color list must equal or longer than t_eval list)
                'b' as blue
                'g' as green
                'r' as red
                'c' as cyan
                'm' as magenta
                'y' as yellow
                'k' as black
    """
    Labels_X1 = 'z-position [\u03BCm]'  # mu: \u03BC
    Labels_X2 = 'Time [sec]'  # mu: \u03BC
    Labels_Y = 'Concentration [mM]'
    Title = ['BDF - Bounded - Diffusion at various time', 'BDF - Unbounded - Diffusion at various time',
             'BDF - Bounded + Unbounded - Diffusion at various time']
    color = ['b-', 'r-', 'g-', 'y-', 'k-', 'm-', 'c-']

    Robin_sim = FVM_1D_Sim_Dslither_Dl_Dr(Ns, Nz, km, K, ka, kd, Dpore_l, Dpore_r, Dslither, mu, rho, Q_r, Qt, Rg, Lp, Rp0,
                                          GAM_0, C0, v0, Bmax, t_array, t_eval, Boundary_Condition, Dirichlet_Conc,
                                          Neumann_Flux, Robin_km, Robin_Conc, Labels_X1, Labels_Y, Title, color)

    """
    Call the result plot from class FVM_1D_Sim:
        output includes plot and solution data from ivp solver
    """
    FVM_Robin_Results = Robin_sim.get_Results()
    Robin_Y1_Data = np.real(FVM_Robin_Results.y)

    model_data = np.zeros([len(t_eval), Nz])
    for j in range(Nz):
        model_data[:, j] = Robin_Y1_Data[j, :]+Robin_Y1_Data[j+Nz, :]

    xx = np.linspace(0, Lp, Nz+1, endpoint=True)
    xarray = np.tile(xx, Ns)

    'Get element center points'
    'Compute the location of the midpoints for each element.'
    xmiddle = (xarray[1:Nz + 1] + xarray[0:Nz])/2
    xmid = np.tile(xmiddle, Ns)
    X1_Data = xmid[:Nz]*10**4

    err = 0
    for i in range(np.shape(Images_data.iloc[:, :])[0]):
        diff = Robin_Y1_Data[:Nz, i]+Robin_Y1_Data[Nz:, i] - \
            Images_data.iloc[i, :]
        err += np.sum(diff**2)
    err_vec.append(err)

    'Plot 1: Result Figure'
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            linestyle = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5))]
            for i in range(len(pt)):
                FVM_Robin_py, = ax.plot(
                    X1_Data, (Robin_Y1_Data[:Nz, pt[i]]+Robin_Y1_Data[Nz:, pt[i]])*1e6, color='black', linestyle=linestyle[i])
                Crystal_1_Images, = ax.plot(
                    X1_Data, Images_data.iloc[pt[i], :]*1e6, color='green', linestyle=linestyle[i], label=f"t={t_eval[pt[i]]} sec")

            line1 = Line2D([0, 1], [0, 1], linestyle='-', color='black')
            line2 = Line2D([0, 1], [0, 1], linestyle='-', color='g')

            plt.xlabel(Labels_X1)
            plt.ylabel(Labels_Y)
            ax.tick_params(direction="in", length=4, width=3)
            lines = plt.gca().get_lines()
            leg1 = plt.legend([lines[i] for i in [1, 3, 5, 7, 9]], [lines[i].get_label() for i in [1, 3, 5, 7, 9]],  loc='lower left', bbox_to_anchor=(
                1, 0.15), shadow=True, mode=None, fancybox=True)
            leg2 = plt.legend([line1, line2], ['model', 'data'],
                              loc='upper left', bbox_to_anchor=(1, 1), shadow=True,
                              mode=None, fancybox=True)
            plt.gca().add_artist(leg1)
            plt.show()
            fig.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_result_all.svg',
                         format='svg', dpi=400,  bbox_inches='tight')

    'Plot 2: data taken just inside the crystal at the left boundary x = 0'
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            Robin_Y1_Data_BC = np.zeros((2, len(t_eval)))
            Images_data_BC = np.zeros((len(t_eval)))

            for i in range(np.shape(Robin_Y1_Data[:, :])[1]):
                Robin_Y1_Data_BC[0, i] = Robin_Y1_Data[0, i]
                Robin_Y1_Data_BC[1, i] = Robin_Y1_Data[Nz, i]
                Images_data_BC[i] = Images_data.iloc[i, 0]
            BCs_data_C_q,      = ax.plot(
                t_eval, (Robin_Y1_Data_BC[0]+Robin_Y1_Data_BC[1])*1e6, 'r--')
            BCs_data_images,   = ax.plot(
                t_eval, Images_data_BC*1e6, 'r-')
            plt.xlabel(Labels_X2)
            plt.ylabel(Labels_Y)
            ax.tick_params(direction="in", length=4, width=3)
            plt.show()
            fig.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_left.svg',
                         format='svg', dpi=400,  bbox_inches='tight')
            
    'Plot 3: data taken just inside the crystal at the right boundary x = Lp'
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            Robin_Y1_Data_BC = np.zeros((2, len(t_eval)))
            Images_data_BC = np.zeros((len(t_eval)))

            for i in range(np.shape(Robin_Y1_Data[:, :])[1]):
                Robin_Y1_Data_BC[0, i] = Robin_Y1_Data[Nz-1, i]
                Robin_Y1_Data_BC[1, i] = Robin_Y1_Data[-1, i]
                Images_data_BC[i] = Images_data.iloc[i, -1]

            BCs_data_C_q,      = ax.plot(
                t_eval, (Robin_Y1_Data_BC[0]+Robin_Y1_Data_BC[1])*1e6, 'r--')
            BCs_data_images,   = ax.plot(
                t_eval, Images_data_BC*1e6, 'r-')

            plt.xlabel(Labels_X2)
            plt.ylabel(Labels_Y)
            ax.tick_params(direction="in", length=4, width=3)
            plt.show()
            fig.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_right.svg',
                         format='svg', dpi=400,  bbox_inches='tight')

    'Plot 4: data taken in the middle of the crystal x = L/2'
    for _ in xtal_fold_to_run:
        if xtal_fold_i == _:
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            Robin_Y1_Data_mid = np.zeros((2, len(t_eval)))
            Images_data_mid = np.zeros((len(t_eval)))

            for i in range(np.shape(Robin_Y1_Data[:, :])[1]):
                Robin_Y1_Data_mid[0, i] = Robin_Y1_Data[int(Nz/2 - 1), i]
                Robin_Y1_Data_mid[1, i] = Robin_Y1_Data[int(Nz + Nz/2 - 1), i]
                Images_data_mid[i] = Images_data.iloc[i, int(Nz/2 - 1)]
            mid_data_C_q,      = ax.plot(
                t_eval, (Robin_Y1_Data_mid[0]+Robin_Y1_Data_mid[1])*1e6, 'r--')
            mid_data_images,   = ax.plot(
                t_eval, Images_data_mid*1e6, 'r-')

            plt.xlabel(Labels_X2)
            plt.ylabel(Labels_Y)
            ax.tick_params(direction="in", length=4, width=3)
            leg = plt.legend([BCs_data_C_q, BCs_data_images], [
                'model', 'data'], loc='upper left', bbox_to_anchor=(
                1, 1), shadow=True, mode=None, fancybox=True)
            plt.show()
            fig.savefig(save_file_path +str(xtal)+'_stripe'+str(xtal_fold_i)+'_middle.svg',
                         format='svg', dpi=400,  bbox_inches='tight')
           
print('SSD = ')
for _ in err_vec:
    print('{:.2e}'.format(_))

# %%
