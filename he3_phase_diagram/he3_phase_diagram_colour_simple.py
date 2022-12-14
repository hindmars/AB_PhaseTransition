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

savefig = True

fill_col_norm = '#B1DFE4'
fill_col_A = '#F5EC20'
fill_col_B = '#CCE3A7'
fill_col_B_super_constP = '#F9F5A9'
fill_col_B_super_varyP = '#CCE3A7'

p_a = (0,40)
T_a = (0, 2.5)

T_smooth = np.linspace(np.min(T_a), np.max(T_a), 200)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)

Tc_line = h.Tc_mK_expt(p_smooth)
TAB_line = h.TAB_mK_expt(p_smooth)

p_melt_line = h.p_melt(T_smooth)

# p_coarse = p_smooth[::10]



#%%
fig, ax = plt.subplots(figsize=(4,3))

p_smooth_sf = p_smooth[p_smooth < h.p_A_bar]
Tc_line_sf = Tc_line[p_smooth < h.p_A_bar]
TAB_line_sf = TAB_line[p_smooth < h.p_A_bar]

# Tc_line_n = Tc_line[p_smooth >= h.p_A_bar]
T_smooth_n = T_smooth[T_smooth>=h.Tc_mK_expt(h.p_A_bar)]



ax.plot(Tc_line_sf, p_smooth_sf, 'k', label=r'$T_c$')
ax.plot(TAB_line_sf, p_smooth_sf, 'k--', label=r'$T_{AB}$')
ax.plot(T_smooth, p_melt_line, 'k-.', label=r'$T_{\rm melt}$')


#%%
# Fillss
p_upper = p_a[1]

#B phase

ax.fill_between(T_a, [h.p_A_bar]*2, color=fill_col_B)



#A phase
boundary_A_lower = np.max(np.stack([Tc_line_sf, TAB_line_sf]),axis=0)
boundary_A_upper = p_upper * np.ones_like(boundary_A_lower)


ax.fill_between(TAB_line_sf, p_smooth_sf, h.p_A_bar * np.ones_like(p_smooth_sf),  color=fill_col_A)
ax.fill_between(Tc_line_sf[p_smooth_sf > h.p_pcp_bar], p_smooth_sf[p_smooth_sf > h.p_pcp_bar], 
                h.p_A_bar * np.ones_like(p_smooth_sf[p_smooth_sf > h.p_pcp_bar]),  color=fill_col_A)

#Normal phase
ax.fill_between(Tc_line_sf, p_a[0] * np.ones_like(p_smooth_sf), p_smooth_sf,  color=fill_col_norm)
ax.fill_between(T_smooth_n, p_a[0] * np.ones_like(T_smooth_n),  h.p_A_bar * np.ones_like(T_smooth_n), color=fill_col_norm)

#%%
ax.set_ylim(*p_a)
ax.set_xlim(*T_a)
ax.grid()
ax.set_xlabel(r'$T$/mK')
ax.set_ylabel(r'$p$/bar')
ax.set_yticks(np.arange(*p_a, 4))
# ax.legend()

# ax.annotate("B phase",(1.8,26), fontsize="larger")
# ax.annotate("A phase",(2.1,29), fontsize="larger")
# ax.annotate("Normal",(2.25,18), fontsize="larger")
# ax.annotate("Solid",(2.0,35), fontsize="larger")
ax.annotate("B phase",(1.0,20), fontsize="larger")
ax.annotate("A",(2.15,29), fontsize="larger")
ax.annotate("Normal",(1.75,4), fontsize="larger")
ax.annotate("Solid",(1.0,36), fontsize="larger")

ax.set_title(r'$^3$He bulk phase diagram, $B=0$')

fig.tight_layout()

if savefig:
    fig.savefig("he3_phase_diagram_simple.pdf")
