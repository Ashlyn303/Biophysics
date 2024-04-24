# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:28:20 2021

@author: Ashlyn
"""
import numpy as np
import math


def FVM_Matrices(D_of_current_element, D_of_previous_element, D_of_next_element, Gam0, Gam1A, Gam1B, Q_x, r, L_of_current_element, L_of_previous_element, L_of_next_element, x_at_left_edge):

    t2 = r**2
    t3 = L_of_next_element + L_of_current_element
    t4 = L_of_current_element + L_of_previous_element
    t5 = L_of_next_element/2.0
    t6 = L_of_current_element/2.0
    t7 = L_of_previous_element/2.0
    t10 = t5+t6
    t11 = t6+t7
    t12 = 1.0/t10
    t13 = 1.0/t11
    Accum_Jac = L_of_current_element*t2*math.pi
    Advec_Jac = np.array([-Q_x*(t7*t13-1.0), Q_x*(t6*t12-1.0) +
                         Q_x*t7*t13, L_of_current_element*Q_x*t12*(-1.0/2.0)])
    Difsn_Jac = np.array([D_of_previous_element*t2*math.pi * (2.0/t4), D_of_current_element *
                         t2*math.pi * -((2.0/t3)+(2.0/t4)), D_of_next_element*t2*math.pi * (2.0/t3)])
    GAM1A_Jac = Gam1A*L_of_current_element*t2*math.pi
    GAM1B_Jac = Gam1B*L_of_current_element*t2*math.pi
    Gam_0_B = Gam0*L_of_current_element*t2*math.pi

    return Accum_Jac, Advec_Jac, Difsn_Jac, GAM1A_Jac, GAM1B_Jac, Gam_0_B

# %%
