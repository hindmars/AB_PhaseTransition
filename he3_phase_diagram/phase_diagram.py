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

h.DEFAULT_T_SCALE = 'PLTS'

p_a = (21,34)
p = np.linspace(*p_a,500,endpoint=True)
p_nodes = range(0,35)

fig, ax = plt.subplots()

Tc_exp = h.Tc_mK_expt(p)
TAB_exp = h.TAB_mK_expt(p)

ax.plot(Tc_exp, p, 'k', label=r"$T_{c}$, Greywall 1986")
ax.plot(TAB_exp, p, 'k--', label=r"$T_{AB}$, Greywall 1986")

def do_TAB_plot(ax, sc_type, p):
    h.DEFAULT_SC_CORRS = sc_type    
    TAB = h.t_AB(p)*h.Tc_mK(p)
    ax.plot(TAB, p, label=r"$T_{AB}$," + sc_type)

sc_type_list = ["Choi-poly", "WS15", "RWS19", "Wiman-thesis"]

for sc_type in sc_type_list:
    do_TAB_plot(ax, sc_type, p)

sc_type = 'Choi-interp'
h.DEFAULT_SC_CORRS = sc_type    
TAB = h.t_AB(p_nodes)*h.Tc_mK(p_nodes)
ax.plot(TAB, p_nodes, 'x', label=r"$T_{AB}$," + sc_type)

ax.grid()
ax.legend(loc='upper right')
ax.set_xlim(1.8, 2.5)
ax.set_ylim(*p_a)
# plt.ylim(0,34)
ax.set_xlabel(r"$T$/mK")
ax.set_ylabel(r"$p$/bar")
ax.set_title(r"$T_{AB}$ predicted vs measured, T scale " + h.DEFAULT_T_SCALE)


