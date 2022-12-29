#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 18:48:46 2022

@author: hindmars
"""

import numpy as np
# import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

import confine as c

rng = np.random.default_rng()

N = 20
R = 10

Areal = rng.standard_normal((N,3,3))
Aimag = rng.standard_normal((N,3,3))

p = 30
t1 = h.t_AB(p)
t2 = 1
t = (t1+t2)/2

eps = 0.2 * h.delta_A_norm(t, p) 
A_init = eps*(Areal + 1j*Aimag)


pot = hw.quartic_potential(t,p)
gr = hw.grid_1d(N, R)

t_eval = np.linspace(0,200,6)
sol_list = hw.relax_from_ic(t_eval, A_init, pot, gr)

#%%


c.plot_confine(*sol_list[1])
c.plot_confine(*sol_list[-1])
