#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 18:21:36 2022

@author: hindmars
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle as pic
import he3_tools as h

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")


bc = 'min'
w = 1100

id = "_29x25"

with open(f'confine_sol_data_{bc:}_1100_BB{id:}.pickle', 'rb') as f:
    Apg_BB_list, s_bw_BB_arr, s_tot_BB_arr, T_a, p_a = pic.load(f)
    
with open(f'confine_sol_data_{bc:}_1100_AB{id:}.pickle', 'rb') as f:
    Apg_AB_list, s_bw_AB_arr, s_tot_AB_arr, T_a, p_a = pic.load(f)
    
# with open(f'confine_sol_data_{bc:}_1100_BA.pickle', 'rb') as f:
#     Apg_BA_list, s_bw_BA_arr, s_tot_BA_arr, T_a, p_a = pic.load(f)
    
with open(f'confine_sol_data_{bc:}_1100_AA{id:}.pickle', 'rb') as f:
    Apg_AA_list, s_bw_AA_arr, s_tot_AA_arr, T_a, p_a = pic.load(f)
    

ds_bw_BB_AB = s_bw_BB_arr[:,:] - s_bw_AB_arr[:,:]
ds_bw_BB_AA = s_bw_BB_arr[:,:] - s_bw_AA_arr[:,:]
ds_tot = s_tot_BB_arr[:,:,0] - s_tot_AA_arr[:,:,0]

#%%

def get_levels_cols(min_level, max_level, n_level):
    cmap_neg = plt.cm.get_cmap("Blues")
    cmap_pos = plt.cm.get_cmap("Reds")
    levels = np.linspace(min_level, max_level, n_level)
    diff_level = (max_level - min_level)/n_level
    cols = list( cmap_neg((levels[levels <0]-diff_level)/(min_level-diff_level)) ) + list(cmap_pos((levels[levels >= 0]+diff_level)/(max_level+diff_level)))
    return levels, cols

#%%

fig, ax = plt.subplots()

levels, cols = get_levels_cols(-1.0, 1.2, 23)

cs = ax.contourf(T_a, p_a, ds_tot, levels, colors=cols)

cbar = fig.colorbar(cs)
cbar.ax.set_ylabel(r'$(\Sigma_{BB}-\Sigma_{AA})/\xi_{GL}(T)|f_B|$')

ax.plot(h.Tc_mK_expt(p_a), p_a, 'k', label=r'$T_c$')
ax.plot(h.TAB_mK_expt(p_a), p_a, 'k--', label=r'$T_{AB}$')

ax.grid()
ax.set_xlabel(r'$T$/mK')
ax.set_ylabel(r'$p$/bar')
ax.legend()
ax.set_title(f'Total surface energy difference, $D = {w:}$ nm \n (wall,bulk) = AA or BB, {bc:}imal pair breaking')

#%%

fig_bw_BB_AA, ax_bw_BB_AA = plt.subplots()

if bc == "min":
    min_l = 0.0
    max_l = 0.9
    n_l = 19
elif bc == "max":
    min_l = 0.0
    max_l = 0.6
    n_l = 13

levels_bw_BB_AA, cols_bw_BB_AA = get_levels_cols(min_l,max_l,n_l)


cs_bw_BB_AA = ax_bw_BB_AA.contourf(T_a, p_a, ds_bw_BB_AA, levels_bw_BB_AA, colors=cols_bw_BB_AA)

cbar_bw_BB_AA = fig_bw_BB_AA.colorbar(cs_bw_BB_AA)
cbar_bw_BB_AA.ax.set_ylabel(r'$(\sigma_{BB} - \sigma_{AA})/\xi_{GL}(T)|f_B|$')

ax_bw_BB_AA.plot(h.Tc_mK_expt(p_a), p_a, 'k', label=r'$T_c$')
ax_bw_BB_AA.plot(h.TAB_mK_expt(p_a), p_a, 'k--', label=r'$T_{AB}$')

ax_bw_BB_AA.grid()
ax_bw_BB_AA.set_xlabel(r'$T$/mK')
ax_bw_BB_AA.set_ylabel(r'$p$/bar')
ax_bw_BB_AA.legend()
ax_bw_BB_AA.set_title(f'Boundary surface energy difference: \n (wall,bulk) = AA or BB, {bc:}imal pair breaking')


#%%

fig_bw_BB_AB, ax_bw_BB_AB = plt.subplots()

if bc == "min":
    min_l = -0.6
    max_l = 0.1
    n_l = 15
elif bc == "max":
    min_l = -3e-3
    max_l = 3e-3
    n_l = 7

levels_bw_BB_AB, cols_bw_BB_AB = get_levels_cols(min_l,max_l,n_l)

cs_bw_BB_AB = ax_bw_BB_AB.contourf(T_a, p_a, ds_bw_BB_AB, levels_bw_BB_AB, colors=cols_bw_BB_AB)

cbar_bw_BB_AB = fig_bw_BB_AB.colorbar(cs_bw_BB_AB)
cbar_bw_BB_AB.ax.set_ylabel(r'$(\sigma_{BB} - \sigma_{AB})/\xi_{GL}(T)|f_B|$')

ax_bw_BB_AB.plot(h.Tc_mK_expt(p_a), p_a, 'k', label=r'$T_c$')
ax_bw_BB_AB.plot(h.TAB_mK_expt(p_a), p_a, 'k--', label=r'$T_{AB}$')

ax_bw_BB_AB.grid()
ax_bw_BB_AB.set_xlabel(r'$T$/mK')
ax_bw_BB_AB.set_ylabel(r'$p$/bar')
ax_bw_BB_AB.legend()
ax_bw_BB_AB.set_title(f'Boundary surface energy difference: \n (wall,bulk) = AB or BB, {bc:}imal pair breaking')
