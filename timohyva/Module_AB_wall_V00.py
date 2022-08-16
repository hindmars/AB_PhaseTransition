#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:23:10 2022

Balanced He3 phase boundary solutions with T = TAB(p)

@author: hindmars, re-edited and modfied by Quang (timohyva@github)

This quisi-newtonian iteration from Scipy lib can't converge when pressure lower than 21 bar.
And it's T-independent, then it will be hard to converge when close to PCP.

"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

# p = 25.5
# p = 25

def get_wall_n_surfaceEnergy(p):

    t = h.t_AB(p);print("\n t_AB is ",t)

    A, pot, gr = hw.krylov_bubble(t,p, gr_pars=(500,125), dim=1)
    # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
    
    # sigma_AB = hw.energy(A,pot,gr)[0]*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))

    # sigma_AB = (hw.energy(*Apg_kry)[0]*h.xi(0,p)/(abs(Apg_kry[1].mat_pars.f_B_norm())*h.xi(t,p)))/np.sqrt(5./3.)

    sigma_AB = (hw.energy(A,pot,gr)[0]*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p)))/np.sqrt(5./3.)

    print('\nSurface energy on event-pressure',p,'is', sigma_AB)
    
    # return (A,pot,gr), sigma_AB
    return sigma_AB



