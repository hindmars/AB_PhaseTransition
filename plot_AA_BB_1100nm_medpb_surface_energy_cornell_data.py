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


bc = 'med'
w = 1100

id = "_29x25"

with open(f'confine_sol_data_{bc:}_1100_BB{id:}.pickle', 'rb') as f:
    Apg_BB_list, s_bw_BB_arr, s_tot_BB_arr, T_a, p_a = pic.load(f)
    
# with open(f'confine_sol_data_{bc:}_1100_AB{id:}.pickle', 'rb') as f:
#     Apg_AB_list, s_bw_AB_arr, s_tot_AB_arr, T_a, p_a = pic.load(f)
    
# with open(f'confine_sol_data_{bc:}_1100_BA.pickle', 'rb') as f:
#     Apg_BA_list, s_bw_BA_arr, s_tot_BA_arr, T_a, p_a = pic.load(f)
    
with open(f'confine_sol_data_{bc:}_1100_AA{id:}.pickle', 'rb') as f:
    Apg_AA_list, s_bw_AA_arr, s_tot_AA_arr, T_a, p_a = pic.load(f)
    

# ds_bw_BB_AB = s_bw_BB_arr[:,:] - s_bw_AB_arr[:,:]
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

# levels, cols = get_levels_cols(-1.0, 1.2, 23)
levels, cols = get_levels_cols(-0.5, 1.0, 16)

cs = ax.contourf(T_a, p_a, ds_tot, levels, colors=cols)

cbar = fig.colorbar(cs)
cbar.ax.set_ylabel(r'$(\Sigma_{BB}-\Sigma_{AA})/\xi_{GL}(T)|f_B|$')

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

ax.grid()
ax.set_xlabel(r'$T$/mK')
ax.set_ylabel(r'$p$/bar')
ax.legend()
ax.set_title(f'Total surface energy difference, $D = {w:}$ nm \n (wall,bulk) = (A,A) or (B,B) {bc:}ium pair breaking')

#%%
fig.savefig(f'surface_energy_AA_BB_1100nm_{bc:}pb{id:}_cornell_data.pdf')
