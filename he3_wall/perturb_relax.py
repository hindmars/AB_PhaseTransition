#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:00:25 2022

For testing stability of solutions found by Newton-Krylov method.

@author: hindmars
"""

import numpy as np
import pickle as pic
import he3_tools as h
import he3_wall as hw
import confine as c

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")


bc = 'max'
w = 1100
id = "_29x25"
id2 = "_Atwist"

with open(f'confine_sol_data_{bc:}_1100_AA{id2:}.pickle', 'rb') as f:
    Apg_AA_list, s_bw_AA_arr, s_tot_AA_arr, T_a, p_a = pic.load(f)

if bc == 'max':
    bcs = [hw.bc_max_pb]*2
elif bc == 'min':
    bcs = [hw.bc_min_pb]*2
    
#%%
m=10
n=10
eps = 0.2

A_init = Apg_AA_list[m][n][0]
pot = Apg_AA_list[m][n][1]
gr = Apg_AA_list[m][n][2]

# Perturb towards Ay and away from Ayy
t_eval = np.linspace(0, 15, 4)

D = h.D_dict["Ay"]- h.D_dict["Ayy"]

perturb = eps * np.multiply.outer(np.sin(np.pi*gr.x/gr.R), D)

Apg_relaxed = hw.relax_from_ic(t_eval, A_init + perturb, pot, gr, bcs)

#%%
ax_list = []

for n, t in enumerate(t_eval):
    ax_list.append(c.plot_confine(*Apg_relaxed[n], bcs))
    ax_list[-1][0].annotate(r"Relaxation time $\tau = {:.1f}$".format(t), (15,2))
    ax_list[-1][0].get_figure().savefig(f'confine_perturb_{bc:}_1100_AA{id2:}_tau{t:.0f}.pdf')

