#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 12:34:57 2022

Plot the total surface energy difference for a channel of 1100 nm, 
along with Cornell group data
(T,p) for AB transition in Heat Exchanger Chamber (HEC) and Isolated Chamber (IC)
- constant pressure coolinng
- varied pressure cooling (reduce from 29.3 bar)

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle as pic
import he3_tools as h

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")


bc = 'max'
w = 1100
id = "_29x25"
id2 = "_Asplay"

with open(f'confine_sol_data_{bc:}_1100_BB{id:}.pickle', 'rb') as f:
    Apg_BB_list, s_bw_BB_arr, s_tot_BB_arr, T_a, p_a = pic.load(f)
    
# with open(f'confine_sol_data_{bc:}_1100_AB{id:}.pickle', 'rb') as f:
#     Apg_AB_list, s_bw_AB_arr, s_tot_AB_arr, T_a, p_a = pic.load(f)
    
# with open(f'confine_sol_data_{bc:}_1100_BA{id:}.pickle', 'rb') as f:
#     Apg_BA_list, s_bw_BA_arr, s_tot_BA_arr, T_a, p_a = pic.load(f)
    
with open(f'confine_sol_data_{bc:}_1100_AA{id2:}.pickle', 'rb') as f:
    Apg_AA_list, s_bw_AA_arr, s_tot_AA_arr, T_a, p_a = pic.load(f)
    

# ds_bw_BB_AB = s_bw_BB_arr[:,:] - s_bw_AB_arr[:,:]
# ds_bw_BB_AA = s_bw_BB_arr[:,:] - s_bw_AA_arr[:,:]

# AA with splay is calculated over the whole channel, so divide by two.
ds_tot = s_tot_BB_arr[:,:,0] - s_tot_AA_arr[:,:,0]/2

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

levels, cols = get_levels_cols(-1.3, 0.2, 16)

cs = ax.contourf(T_a, p_a, ds_tot, levels, colors=cols)

cbar = fig.colorbar(cs)
cbar.ax.set_ylabel(r'$(\Sigma_{AA}-\Sigma_{BB})/\xi_{GL}(T)|f_B|$')

p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)
ax.plot(h.Tc_mK_expt(p_smooth), p_smooth, 'k', label=r'$T_c$')
ax.plot(h.TAB_mK_expt(p_smooth), p_smooth, 'k--', label=r'$T_{AB}$')

#%%

# Plot Cornell grpooup data 2022
# first constant pressure

data_constp = np.genfromtxt('Cornell22_const_P.csv', skip_header=1, delimiter=',')
p_data_constp = data_constp[:,0]
T_IC_data_constp = data_constp[:,1]
T_HEC_data_constp = data_constp[:,2]

ax.scatter(T_IC_data_constp, p_data_constp, marker='<', s=15, c='xkcd:lightblue', alpha=0.5, edgecolors='b', linewidth=1)
ax.scatter(T_HEC_data_constp, p_data_constp, marker='<', s=20, c='xkcd:pink', alpha=0.5, edgecolors='k', linewidth=1)

#%%

# Plot Cornell group data 2022
# now varied pressure

data_variedp = np.genfromtxt('Cornell22_varied_P.csv', skip_header=1, delimiter=',')
p_data_variedp = data_variedp[:,0]
T_HEC_data_variedp = data_variedp[:,1]
T_IC_data_variedp = data_variedp[:,2]

ax.scatter(T_IC_data_variedp, p_data_variedp, marker='x', s=15, c='b', linewidth=1)
ax.scatter(T_HEC_data_variedp, p_data_variedp, marker='s', s=15, c='xkcd:pink', alpha=0.5, edgecolors='k', linewidth=1)

#%%

# # Plot Alberta group 2020 confinment transitionn data

# data_shook_et_al = np.genfromtxt('Shook_et_al_20_Fig4_wpd.csv', skip_header=1, delimiter=',')
# t_filled_shook = data_shook_et_al[:,0]
# p_filled_shook = data_shook_et_al[:,1]
# t_empty_shook = data_shook_et_al[:,2]
# p_empty_shook = data_shook_et_al[:,3]

# ax.scatter(t_filled_shook, p_filled_shook, marker='o', s=15, c='xkcd:green', 
#            linewidth=1, label='Shook + 2020')
# ax.scatter(t_empty_shook, p_empty_shook, marker='o', s=15, c='xkcd:white', 
#            edgecolors='xkcd:green', linewidth=1, label='Shook et al 2020 (empty)')

#%%
ax.set_xlim(1.7, 2.4)
ax.set_ylim(18, 30)
ax.grid()
ax.set_xlabel(r'$T$/mK')
ax.set_ylabel(r'$p$/bar')
ax.legend(fontsize='smaller')
ax.set_title(f'Total surface energy difference, $D = {w:}$ nm \n (wall,bulk) = (A,A+splay) or (B,B) {bc:}imal pair breaking', fontsize='smaller')

#%%
fig.savefig(f'surface_energy_AA_BB_1100nm_{bc:}pb{id2:}_cornell_data.pdf')
