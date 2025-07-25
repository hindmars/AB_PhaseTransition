#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:23:10 2022

He3 phase boundary solutions

@author: hindmars
"""

import numpy as np
import he3_tools as h
import he3_wall as hw

h.set_default("DEFAULT_T_SCALE", "Greywall")

#%% Set pressure and rduced temperature

# p = 25.5
p = 22
# t = 1.5/h.Tc_mK(p)
t = h.tAB(p)

#%% Lattice size: grid points and physical llength

N, L_xi = (500, 25)
savefig=False

#%% Get wall solution

Apg_kry, sigma_bw, sigma_tot = hw.get_wall(t, p, L_xi * h.xi(t,p), right_phase='B')

# Apg_kry is a tuple
# Apg_kry[0] is order parameter array, Apg_kry[0].shape = (N, 3, 3), complex
# Apg_kry[1] is a quartic potential object
# Apg_kry[2] is a grid obkect

#%% Plot wall

ax2 = hw.plot_wall(*Apg_kry, plot_gap=True, phase_marker=True)

#%% Plot R1

x = Apg_kry[2].x
R_vec = h.R_terms(Apg_kry[0])
xiGL0 = h.xi(0,p)

ax2[2].plot((x - np.max(x)/2)/xiGL0, 4*R_vec[0]*(1 - R_vec[0]), '--', label=r'$4R_1(1-R_1)$')
ax2[2].legend(fontsize='smaller', bbox_to_anchor=(0.9, 0.5))

#%% Make a figure if wanted

if savefig:
    file_name = 'AB_wall_t={:.2f}_p={:.1f}.pdf'.format(t,p)
    print('print to figure',file_name)
    ax[0].get_figure().savefig(file_name)

