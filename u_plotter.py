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

v = np.linspace(0,1,100)*np.sqrt(2)

delta_A = h.delta_A_norm(t,p)
delta_B = h.delta_B_norm(t,p)
delta_planar = h.delta_planar_norm(t,p)

A_A = np.multiply.outer( v , h.D_A * delta_A)
A_B = np.multiply.outer( v , h.D_B * delta_B)
A_AB = np.multiply.outer( (1 - v) , h.D_A * delta_A) + np.multiply.outer( v , h.D_B * delta_B)
A_Aplanar = np.multiply.outer( (1 - v) , h.D_A * delta_A) + np.multiply.outer( v , h.D_planar * delta_planar)

U_A = h.U(A_A, (t-1), h.beta_norm_asarray(t, p) )
U_B = h.U(A_B, (t-1), h.beta_norm_asarray(t, p) )
U_AB = h.U(A_AB, (t-1), h.beta_norm_asarray(t, p) )
U_Aplanar = h.U(A_Aplanar, (t-1), h.beta_norm_asarray(t, p) )

# plt.plot(v,  U_A)
# plt.plot(v,  U_B)
plt.plot(v,  U_AB, label='X = B')
plt.plot(v,  U_Aplanar, label='X = planar')

plt.xlabel(r'$(\Delta - \Delta_A)/(\Delta_X - \Delta_A)$')
plt.ylabel(r'$\tilde{f}$')
plt.xlim(0,1.4)
plt.legend()
plt.grid()
