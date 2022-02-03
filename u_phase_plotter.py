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

plt.figure()

for phase in h.inert_phases:
    delta = h.delta_phase_norm(t, p, phase)
    D = h.D_dict[phase]
    A = np.multiply.outer( v , D * delta)
    UA = h.U(A, (t-1), h.beta_norm_asarray(t, p) )
    plt.plot(v,  UA, label=phase)

plt.xlabel(r'$\Delta /\Delta_X$')
plt.ylabel(r'$f / [\frac{1}{3} N(0) (k_BT_c)^{2} ]$ ($T/T_c = 0.8$, $p = 26$ bar)')
plt.xlim(0, 1.6)
plt.ylim(-1,0.5)
plt.legend(title="Phase X", loc="lower left")
plt.grid()
