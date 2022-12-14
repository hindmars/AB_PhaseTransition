#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 19:02:45 2022

Check that adjustments to SC coupling corrections really do put the theoretical 
TAB line (RQS19) on top of the experimental one (Greywall)


@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

# h.DEFAULT_T_SCALE = 'PLTS'

savefig = True
p_a = (20,34)
p = np.linspace(*p_a,500,endpoint=True)
p_nodes = range(0,35)

sc_type_list = ["Choi-poly", "WS15", "RWS19", "Wiman-thesis"]

symbol_list = ['-.', '--', '-', ':']


fig, ax = plt.subplots()

Tc_exp = h.Tc_mK_expt(p)
TAB_exp = h.TAB_mK_expt(p)

ax.plot(Tc_exp, p, 'b', label=r"$T_{c}$, Greywall 1986")
ax.plot(TAB_exp, p, 'g', label=r"$T_{AB}$, Greywall 1986")

def do_TAB_plot(ax, sc_type, p, ls_mk='-'):
    h.DEFAULT_SC_CORRS = sc_type    
    TAB = h.t_AB(p)*h.Tc_mK(p)
    ax.plot(TAB, p, ls_mk, c='k', label=r"$T_{AB}$," + sc_type)

sc_type = 'Choi-interp'
h.DEFAULT_SC_CORRS = sc_type    
TAB = h.t_AB(p_nodes)*h.Tc_mK(p_nodes)
ax.plot(TAB, p_nodes, 'x', c='k', label=r"$T_{AB}$," + sc_type)

for sc_type, symb in zip(sc_type_list, symbol_list):
    do_TAB_plot(ax, sc_type, p, symb)

ax.grid()
ax.legend(loc='upper right')
ax.set_xlim(1.8, 2.5)
ax.set_ylim(*p_a)
# plt.ylim(0,34)
ax.set_xlabel(r"$T$/mK")
ax.set_ylabel(r"$p$/bar")
ax.set_title(r"$T_{AB}$ predicted vs measured, T scale " + h.DEFAULT_T_SCALE)


if savefig:
    fig.savefig('he3_tAB_compare_{}.pdf'.format(h.DEFAULT_T_SCALE))