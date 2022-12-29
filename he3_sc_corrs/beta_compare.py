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

savefig = False

cols = ['r', 'g', 'b', 'c', 'm']

lbls = [ r'$\beta_{{ {} }}$'.format(n) for n in range(1,6) ]
lbls_AB = [ r'$\beta_{{ {} }}$'.format(l) for l in ['A', 'B'] ]

p_a = (0,34)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)
p_nodes = np.linspace(0,34,35)

def beta_plot(p, beta_list, ls_mk='-', col_list=cols, leg_list=[None]*5):
    for n, (beta, col) in enumerate(zip(beta_list, col_list)):
        if n+1 in [2,3,4] or len(beta_list) == 2:
            # ax[0].plot(p, beta, ls_mk, c=col, label=r'$\beta_{{ {} }}$'.format(n + 1))
            ax[0].plot(p, beta, ls_mk, c=col, label=leg_list[n])
        elif n+1 in [1,5]:
            ax[1].plot(p, beta, ls_mk, c=col, label=leg_list[n])
    
fig, ax = plt.subplots(2,1, sharex=True)

beta_RWS19_arr = h.beta_norm_asarray(1, p_smooth)
beta_plot(p_smooth, beta_RWS19_arr, '-', cols, lbls)
    
h.set_default('DEFAULT_SC_CORRS', 'WS15')
beta_WS15_arr = h.beta_norm_asarray(1, p_smooth)
beta_plot(p_smooth, beta_WS15_arr, '--', cols)

h.set_default('DEFAULT_SC_CORRS', 'Choi-interp')
beta_Choi_arr = h.beta_norm_asarray(1, p_nodes)
beta_plot(p_nodes, beta_Choi_arr, 'x', cols)

h.set_default('DEFAULT_SC_CORRS', 'Choi-poly')
beta_ChoiPoly_arr = h.beta_norm_asarray(1, p_smooth)
beta_plot(p_smooth, beta_ChoiPoly_arr, '-.', cols)

# for n, (beta, col) in enumerate(zip(beta_ChoiPoly_arr, cols)):
#     plt.plot(p_smooth, np.array(beta), '-', c=col)


for a in ax:
    a.set_xlim(*p_a)
    # a.set_xlim(p_a[0], p_a[1]+10)
    a.grid()
    a.legend(bbox_to_anchor = (1.0, 0.75))

ax[0].set_title(r'$\beta_a$: - RWS19; -- WS15; x Choi; -.; Choi poly')

ax[0].set_ylim(0.01, 0.025)
ax[1].set_ylim(-0.03, -0.005)

ax[1].set_xlabel(r'$p$/bar')
fig.tight_layout()

if savefig:
    fig.savefig('beta_compare.pdf')
