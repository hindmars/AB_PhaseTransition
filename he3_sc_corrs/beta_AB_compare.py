#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 18:02:02 2022

Stron coupling parameters comparator

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import he3_tools as h

savefig = True

# cols = ['r', 'g', 'b', 'c', 'm']
prop_cycle = plt.rcParams['axes.prop_cycle']
cols = prop_cycle.by_key()['color']

lbls = [ r'$\beta_{{ {} }}$'.format(n) for n in range(1,6) ]
lbls_AB = [ r'$\beta^{{\rm sc}}_{{ {} }}/\beta_0$'.format(l) for l in ['A', 'B'] ]

sc_type_list = ["Choi-poly", "WS15", "RWS19", "Wiman-thesis"]

symbol_list = ['-.', '--', '-', ':']


p_a = (0,34)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)
p_nodes = np.linspace(0,34,35)

def beta_plot(ax, p, beta_list, ls_mk='-', col_list=cols, leg_list=[None]*5):
    for n, (beta, col) in enumerate(zip(beta_list, col_list)):
            ax.plot(p, beta, ls_mk, c=col, label=leg_list[n])
    
fig, ax = plt.subplots()

def do_one_plot(ax, sc_type, p, ls_mk, lbls=[None]*5):
    h.set_default('DEFAULT_SC_CORRS', sc_type)
    beta_A = h.beta_A_norm(1, p)/h.beta_const
    beta_B = h.beta_B_norm(1, p)/h.beta_const
    beta_plot(ax, p, [beta_A, beta_B], ls_mk, cols[0:2], lbls)
    

def do_plots(ax):

    do_one_plot(ax, 'Choi-interp', p_nodes, 'x', lbls_AB)

    for sc_type, symbol in zip(sc_type_list, symbol_list):
        do_one_plot(ax, sc_type, p_smooth, symbol)

    ax.axvline(h.p_pcp_bar, ls='--', label=r'$p_{\rm PCP}$')

    return

do_plots(ax)

ax.set_xlim(p_a[0], p_a[1])
# ax.set_ylim(0.011, 0.021)
ax.grid()
ax.legend(loc='lower left')
ax.set_xlabel(r'$p$/bar')
ax.set_title(r'$\beta^{\rm sc}/\beta_0$: x Choi; -. Choi poly; - RWS19; -- WS15; : Wiman thesis')

# Zoom
x1 = 20.5
x2 = 22.5
y1 = 1.25
y2 = 1.35

# Make the zoom-in plot:
axins = zoomed_inset_axes(ax, 5, loc='upper right') # zoom = 2
# axins.plot(p_smooth, beta_ChoiPoly_A)

do_plots(axins)
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# plt.xticks(visible=False)
# plt.yticks(visible=False)
axins.grid()
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

fig.tight_layout()


if savefig:
    fig.savefig('beta_AB_compare.pdf')

