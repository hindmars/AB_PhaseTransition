#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:34:55 2022

@author: hindmars
"""

import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import he3_tools as h

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")

#%%

# marker_data_IC_constP  = {"marker":'<', "s":40, "c":"#0000FF", "alpha":1.0, "edgecolors":'b', "linewidth":1, "label":"IC"}
# marker_data_HEC_constP = {"marker":'<', "s":25, "c":"#FF80FF", "alpha":1.0, "edgecolors":'k', "linewidth":1, "label":"HEC"}

# marker_data_IC_varyP  = {"marker":'x', "s":40, "c":"#0000FF", "alpha":1.0, "edgecolors":'b', "linewidth":1, "label":"IC"}
# marker_data_HEC_varyP = {"marker":'s', "s":40, "c":"#FF80FF", "alpha":1.0, "edgecolors":'k', "linewidth":1, "label":"HEC"}


fill_col_norm = '#B1DFE4'
fill_col_A = '#F5EC20'
fill_col_B = '#CCE3A7'
fill_col_B_super_constP = '#F9F5A9'
fill_col_B_super_varyP = '#CCE3A7'

p_a = (18,30)
T_a = (1.7, 2.4)

T_smooth = np.linspace(np.min(T_a), np.max(T_a), 200)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)

Tc_line = h.Tc_mK_expt(p_smooth)
TAB_line = h.TAB_mK_expt(p_smooth)

h.DEFAULT_SC_CORRS = "RWS19"
h.DEFAULT_SC_ADJUST=False
TAB_RWS = h.t_AB(p_smooth)*h.Tc_mK(p_smooth)

p_coarse = p_smooth[::10]
# TAB_RWS_adj = h.t_AB(p_coarse, sc_adjust=True)*h.Tc_mK(p_coarse)
h.DEFAULT_SC_ADJUST=True
TAB_RWS_adj = h.t_AB(p_coarse)*h.Tc_mK(p_coarse)


#%%
fig, ax = plt.subplots(figsize=(5,3))


ax.plot(Tc_line, p_smooth, 'k', label=r'$T_c$')
ax.plot(TAB_line, p_smooth, 'k', label=r'$T_{AB}$')
ax.plot(TAB_RWS, p_smooth, label=r"$T_{AB}$, RWS 2019")
ax.plot(TAB_RWS_adj, p_coarse, 'x', label=r"$T_{AB}$, RWS 2019, adjusted")


#%%
# Fillss
p_upper = p_a[1]

#B phase

ax.fill_between(T_a, [p_upper]*2, color=fill_col_B)



#A phase
boundary_A_lower = np.max(np.stack([Tc_line, TAB_line]),axis=0)
boundary_A_upper = p_upper * np.ones_like(boundary_A_lower)

ax.fill_between(TAB_line, p_smooth, p_upper * np.ones_like(p_smooth),  color=fill_col_A)
ax.fill_between(Tc_line[p_smooth > h.p_pcp_bar], p_smooth[p_smooth > h.p_pcp_bar], p_upper * np.ones_like(p_smooth[p_smooth > h.p_pcp_bar]),  color=fill_col_A)

#Normal phase
ax.fill_between(Tc_line, p_a[0] * np.ones_like(p_smooth), p_smooth,  color=fill_col_norm)

#%%
ax.set_ylim(*p_a)
ax.set_xlim(1.7, 2.4)
ax.grid()
ax.set_xlabel(r'$T$/mK')
ax.set_ylabel(r'$p$/bar')
ax.legend()

ax.annotate("B phase",(1.85,26), fontsize="larger")
ax.annotate("A phase",(2.15,28), fontsize="larger")
ax.annotate("Normal",(2.25,19), fontsize="larger")

fig.savefig("he3_phase_diagram.pdf")
