#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 18:02:02 2022

Stron coupling parameters comparator

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

savefig = True

cols = ['r', 'g', 'b', 'c', 'm']

lbls = [ r'$\beta_{{ {} }}$'.format(n) for n in range(1,6) ]
lbls_AB = [ r'$\beta_{{ {} }}$'.format(l) for l in ['A', 'B'] ]

p_a = (0,34)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)
p_nodes = np.linspace(0,34,35)

def dbeta_plot(p, dbeta_list, ls_mk='-', col_list=cols, leg_list=[None]*5):
    for n, (dbeta, col) in enumerate(zip(dbeta_list, col_list)):
        ax.plot(p, dbeta, ls_mk, c=col, label=leg_list[n])
    
fig, ax = plt.subplots()
    
dbeta_RWS19_arr = h.delta_beta_norm_asarray(p_smooth)
dbeta_plot(p_smooth, dbeta_RWS19_arr, '-', cols, lbls)
    
h.set_default('DEFAULT_SC_CORRS', 'WS15')
dbeta_WS15_arr = h.delta_beta_norm_asarray(p_smooth)
dbeta_plot(p_smooth, dbeta_WS15_arr, '--', cols)

h.set_default('DEFAULT_SC_CORRS', 'Choi-interp')
dbeta_Choi_arr = h.delta_beta_norm_asarray(p_nodes)
dbeta_plot(p_nodes, dbeta_Choi_arr, 'x', cols)

h.set_default('DEFAULT_SC_CORRS', 'Choi-poly')
dbeta_ChoiPoly_arr = h.delta_beta_norm_asarray(p_smooth)
dbeta_plot(p_smooth, dbeta_ChoiPoly_arr, '-.', cols)

h.set_default('DEFAULT_SC_CORRS', 'Wiman-thesis')
dbeta_WimanThesis_arr = h.delta_beta_norm_asarray(p_smooth)
dbeta_plot(p_smooth, dbeta_WimanThesis_arr, ':', cols)



ax.set_xlim(*p_a)
# a.set_xlim(p_a[0], p_a[1]+10)
ax.grid()
ax.legend(bbox_to_anchor = (1.0, 0.75))

ax.set_title(r'$\Delta\beta_a/\beta_0^{\rm wc}$: x Choi; -. Choi poly; - RWS19; -- WS15;  : Wiman thesis',
             fontsize=10)

ax.set_ylim(-0.8, 0.25)
# ax[1].set_ylim(-0.03, -0.005)

ax.set_xlabel(r'$p$/bar')
fig.tight_layout()

if savefig:
    fig.savefig('sc_corrs_compare.pdf')
