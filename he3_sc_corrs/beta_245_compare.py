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
lbls_345 = r'$\beta^{{\rm sc}}_{{ {} }}/\beta_0$'.format('345') 

sc_type_list = ["Choi-poly", "WS15", "RWS19", "Wiman-thesis"]

symbol_list = ['-.', '--', '-', ':']


p_a = (0,34)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)
p_nodes = np.linspace(0,34,35)

def beta_plot(ax, p, beta, ls_mk='-', col='k', leg=None):
    ax.plot(p, beta, ls_mk, c=col, label=leg)
    return
    
fig, ax = plt.subplots()

def do_one_plot(ax, sc_type, p, ls_mk, lbls=None):
    h.set_default('DEFAULT_SC_CORRS', sc_type)
    beta_245 = (h.beta_norm(1, p, 2) + h.beta_norm(1, p, 4)+ h.beta_norm(1, p, 5)) /h.beta_const
    beta_plot(ax, p, beta_245, ls_mk, 'k', sc_type)
    

def do_plots(ax):

    do_one_plot(ax, 'Choi-interp', p_nodes, 'x')

    for sc_type, symbol in zip(sc_type_list, symbol_list):
        do_one_plot(ax, sc_type, p_smooth, symbol)

    # ax.axvline(h.p_pcp_bar, ls='--', label=r'$p_{PCP}$')

    return

do_plots(ax)

ax.set_xlim(p_a[0], p_a[1])
# ax.set_ylim(0.011, 0.021)
ax.grid()
ax.legend(loc='upper right')
ax.set_xlabel(r'$p$/bar')
ax.set_title(r'$\beta^{\rm sc}_{245}/\beta_0$')

fig.tight_layout()


if savefig:
    fig.savefig('beta_245_compare.pdf')

