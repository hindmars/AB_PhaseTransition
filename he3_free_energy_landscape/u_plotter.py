#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:25:44 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


t = 0.8
p = 26

v = np.linspace(0,1,100)*2

delta_A = h.delta_A_norm(t,p)
delta_B = h.delta_B_norm(t,p)
delta_planar = h.delta_planar_norm(t,p)
delta_polar = h.delta_polar_norm(t,p)

D_AX = np.matmul(h.R_xz, h.D_A)

A_A = np.multiply.outer( v , D_AX * delta_A)
A_B = np.multiply.outer( v , h.D_B * delta_B)
A_AB = np.multiply.outer( (1 - v) , h.D_A * delta_A) + np.multiply.outer( v , h.D_B * delta_B)
A_planar = np.multiply.outer( v , h.D_dict["planar"]*delta_planar)
A_polar = np.multiply.outer( v , h.D_dict["polar"]*delta_polar)

U_A = h.U(A_A, (t-1), h.beta_norm_asarray(t, p) )
U_B = h.U(A_B, (t-1), h.beta_norm_asarray(t, p) )
U_planar = h.U(A_planar, (t-1), h.beta_norm_asarray(t, p) )
U_polar = h.U(A_polar, (t-1), h.beta_norm_asarray(t, p) )


U_AB = h.U(A_AB, (t-1), h.beta_norm_asarray(t, p) )



plt.plot(v,  U_A, label='X = A, Y = N')
plt.plot(v,  U_B, label='X = B, Y = N')
plt.plot(v,  U_planar, label='X = planar, Y = N')
plt.plot(v,  U_polar, label='X = polar, Y = N')
# plt.plot(v,  U_AB, label='X = B, Y = A')


plt.xlabel(r'$(\Delta - \Delta_Y)/(\Delta_X - \Delta_Y)$')
plt.ylabel(r'$f / [\frac{1}{3} N(0) (k_BT_c)^{2} ]$ ($T/T_c = 0.8$, $p = 26$ bar)')
plt.xlim(0, 1.4)
plt.ylim(-1,0.5)
plt.legend()
plt.grid()
