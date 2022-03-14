#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 13:35:33 2021

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


mint, maxt = (0., 1.0)

t = np.linspace(mint, maxt, 100)

ax_list = []

for nax in range(0,2):
    fig, ax = plt.subplots()
    ax_list.append(ax)


vol=1e6**3 # ie mm3

unit_factor = 1e12
unit_str = 'p'


for p in range(20, 34, 2):

    ax_list[0].plot(t*h.T_mK(1,p), (h.f_A_norm(t, p) - h.f_B_norm(t, p)) * h.f_scale(p) * vol * unit_factor, 
                    color='b', alpha=(p-11)/(24), label='{}'.format(p) )
    ax_list[1].plot(t*h.T_mK(1,p), (h.f_A_norm(t, p) - h.f_B_norm(t, p))/np.abs(h.f_B_norm(t, p)), 
                    color='b', alpha=(p-11)/(24), label='{}'.format(p) )




# ax_list[0].xlabel(r'$T/T_c$')
ax_list[0].set_xlabel(r'$T$ / mK')
ax_list[1].set_xlabel(r'$T$ / mK')

# plt.ylabel(r'$\Delta f_{AB} /N(0) (k_B T_c)^2$')
ax_list[0].set_ylabel(r'$\Delta f_{AB}$/ pJ mm$^{-3}$')
ax_list[1].set_ylabel(r'$\Delta f_{AB} / |f_B|$')


# plt.xlim(mint,maxt+0.1)
ax_list[0].set_xlim(1.5,2.5)
ax_list[1].set_xlim(1.5,2.5)

ax_list[0].set_ylim(-5e-4, 5e-4)
ax_list[1].set_ylim(-0.05, 0.05)

ax_list[0].grid()
ax_list[1].grid()

ax_list[0].legend(title=r'$p$/bar')
ax_list[1].legend(title=r'$p$/bar')

ax_list[0].get_figure().tight_layout()
ax_list[1].get_figure().tight_layout()

