#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 19:19:48 2022

Plotting Rcrit_min against R_trans for constant P data Cornell

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")

#%%

# marker_data_IC_constP  = {"marker":'<', "s":40, "c":"#0000FF", "alpha":1.0, "edgecolors":'b', "linewidth":1, "label":"IC"}
# marker_data_HEC_constP = {"marker":'<', "s":25, "c":"#FF80FF", "alpha":1.0, "edgecolors":'k', "linewidth":1, "label":"HEC"}

# marker_data_IC_varyP  = {"marker":'x', "s":40, "c":"#0000FF", "alpha":1.0, "edgecolors":'b', "linewidth":1, "label":"IC"}
# marker_data_HEC_varyP = {"marker":'s', "s":40, "c":"#FF80FF", "alpha":1.0, "edgecolors":'k', "linewidth":1, "label":"HEC"}

# fill_col_norm = '#B1DFE4'
# fill_col_A = '#F5EC20'
# fill_col_B_super_constP = '#F9F5A9'
# fill_col_B_super_varyP = '#CCE3A7'


from parpia_graph_settings import *


#%%
data_constp = np.genfromtxt('Cornell22_const_P.csv', skip_header=1, delimiter=',')

p_constp = data_constp[:,0]
T_IC_constp = data_constp[:,1]
T_HEC_constp = data_constp[:,2]

t_HEC_constp = T_HEC_constp[p_constp > h.p_pcp_bar]/h.Tc_mK(p_constp[p_constp > h.p_pcp_bar])
t_IC_constp = T_IC_constp[p_constp > h.p_pcp_bar]/h.Tc_mK(p_constp[p_constp > h.p_pcp_bar])

Rcrit_HEC = h.critical_radius(t_HEC_constp, p_constp[p_constp > h.p_pcp_bar])
Rcrit_IC = h.critical_radius(t_IC_constp, p_constp[p_constp > h.p_pcp_bar])

Rcrit_min =  []
xi_GL_min = []

N=51

for p in p_constp[p_constp > h.p_pcp_bar]:
    t_path = np.linspace(h.t_AB(p), 1.0 - 1/N, N)
    Rcrit_path = h.critical_radius(t_path, p)
    Rcrit_min.append(np.min(Rcrit_path))
    xi_GL_path = h.xi(t_path, p)
    xi_GL_min.append(np.min(xi_GL_path))

Rcrit_min = np.array(Rcrit_min)
xi_GL_min = np.array(xi_GL_min)

xi_GL_HEC = h.xi(t_HEC_constp, p_constp[p_constp > h.p_pcp_bar])
xi_GL_IC = h.xi(t_IC_constp, p_constp[p_constp > h.p_pcp_bar])


plt.scatter(Rcrit_min, Rcrit_HEC, **marker_data_HEC_constP)
plt.scatter(Rcrit_min, Rcrit_IC, **marker_data_IC_constP)

plt.xlim(0,4e4)
plt.ylim(0,3e4)

plt.xlabel("Minimum AB interface curvature radius, A phase [nm]")
plt.ylabel(r"AB interface curvature radius at transition [nm]")
plt.grid()
plt.legend(loc='upper left')

#55
plt.figure()

plt.scatter(1000/Rcrit_min, 1000/Rcrit_HEC, **marker_data_HEC_constP)
plt.scatter(1000/Rcrit_min, 1000/Rcrit_IC, **marker_data_IC_constP)

plt.xlim(0,0.7)
plt.ylim(0,1.4)
plt.grid()
plt.xlabel(r"Minimum AB interface curvature, A phase [$\mu$m$^{-1}$]")
plt.ylabel(r"AB interface curvature at transition [$\mu$m$^{-1}$]")
plt.legend(loc='upper left')

#%%

plt.figure()

plt.scatter(xi_GL_min/Rcrit_min, xi_GL_HEC/Rcrit_HEC, **marker_data_HEC_varyP)
plt.scatter(xi_GL_min/Rcrit_min, xi_GL_IC/Rcrit_IC, **marker_data_IC_varyP)

plt.xlim(0,0.02)
plt.ylim(0,0.03)
plt.grid()
plt.xlabel(r"$\kappa_{\rm min}^A \xi_{\rm GL}$")
plt.ylabel(r"$\kappa^B \xi_{\rm GL}$")
plt.legend(loc='upper left')
