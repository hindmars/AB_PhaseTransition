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
# p = 25.5
p = 22
# t = 1.5/h.Tc_mK(p)
t = h.tAB(p)

N, L_xi = (500, 25)

savefig=False

# def get_and_plot_wall(t,p, w, **kwargs):

#     (A, pot, gr), sigma_bw, sigma_tot = hw.get_wall(t, p, w, **kwargs)
#     # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
    
#     # Plotting
#     ax = hw.plot_wall(A, pot, gr, plot_gap=True, **kwargs)

#     return (A,pot,gr), ax

#%%

Apg_kry, sigma_bw, sigma_tot = hw.get_wall(t, p, L_xi * h.xi(t,p), right_phase='B')

#%%
# axwall_eigs = hw.plot_eigs(*Apg_kry)

#%% Replot wall

x = Apg_kry[2].x
R_vec = h.R_terms(Apg_kry[0])
xiGL0 = h.xi(0,p)

ax2 = hw.plot_wall(*Apg_kry, plot_gap=True, phase_marker=True)
ax2[2].plot((x - np.max(x)/2)/xiGL0, 4*R_vec[0]*(1 - R_vec[0]), '--', label=r'$4R_1(1-R_1)$')

ax2[2].legend(fontsize='smaller', bbox_to_anchor=(0.9, 0.5))

#%%
# sigma_AB = hw.energy(*Apg_kry)[0]*h.xi(0,p)/(abs(Apg_kry[1].mat_pars.f_B_norm())*h.xi(t,p))

if savefig:
    file_name = 'AB_wall_t={:.2f}_p={:.1f}.pdf'.format(t,p)
    print('print to figure',file_name)
    ax[0].get_figure().savefig(file_name)

