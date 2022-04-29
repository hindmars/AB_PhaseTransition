#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 20:07:00 2022

@author: hindmars
"""

import numpy as np
import pickle as pic
import he3_tools as h
import he3_wall as hw
import confine as c

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")

w = 1100

id = '_Atwist'

T_arr = np.linspace(1.7, 2.4, 29)
p_arr = np.linspace(18, 30, 25)

n_T = len(T_arr)
n_p = len(p_arr)

Apg_conf_BB_list = [ [None for T in T_arr] for p in p_arr ]
Apg_conf_AA_list = [ [None for T in T_arr] for p in p_arr ]

sigma_bw_AA_arr = np.ones((n_p, n_T))*np.nan

sigma_tot_AA_arr = np.ones((n_p, n_T, 3))*np.nan

for m, p in enumerate(p_arr):
    t_arr = T_arr/h.Tc_mK_expt(p)
    for n, t in enumerate(t_arr[t_arr < 1]):

        Apg_conf_AA_list[m][n], sigma_bw_AA_arr[m,n], sigma_tot_AA_arr[m,n] = c.get_confine(t, p, w, "Ay", "Ay", 
                                                    bcs = [hw.bc_max_pb, hw.bc_max_pb], T_list=[h.z3, h.T_xy/2, h.z3], 
                                                    verbose=True, maxiter=50, f_tol=1e-3)
        T = T_arr[n]

        print('AA: P:', p, 'T:', T, ', Surface energy:',  sigma_bw_AA_arr[m,n], 'Total energy:',  sigma_tot_AA_arr[m,n][0])



for m, p in enumerate(p_arr):
    t_arr = T_arr/h.Tc_mK_expt(p)
    for n, t in enumerate(t_arr[t_arr < 1]):
        if Apg_conf_AA_list[m][n] is not None:
            print("after return", m, n, Apg_conf_AA_list[m][n][1].args)

#%%

# confine_sol_data_max_1100_BB = (Apg_conf_BB_list, sigma_bw_BB_arr, sigma_tot_BB_arr, T_arr, p_arr)
# confine_sol_data_max_1100_AB = (Apg_conf_AB_list, sigma_bw_AB_arr, sigma_tot_AB_arr, T_arr, p_arr)
# confine_sol_data_max_1100_BA = (Apg_conf_BA_list, sigma_bw_BA_arr, sigma_tot_BA_arr, T_arr, p_arr)
confine_sol_data_max_1100_AA = (Apg_conf_AA_list, sigma_bw_AA_arr, sigma_tot_AA_arr, T_arr, p_arr)


# with open(f'confine_sol_data_max_1100_BB{id:}.pickle', 'wb') as f:
#     pic.dump(confine_sol_data_max_1100_BB, f)
    
# with open(f'confine_sol_data_max_1100_AB{id:}.pickle', 'wb') as f:
#     pic.dump(confine_sol_data_max_1100_AB, f)
    
# with open(f'confine_sol_data_max_1100_BA{id:}.pickle', 'wb') as f:
#     pic.dump(confine_sol_data_max_1100_BA, f)
    
with open(f'confine_sol_data_max_1100_AA{id:}.pickle', 'wb') as f:
    pic.dump(confine_sol_data_max_1100_AA, f)
    

