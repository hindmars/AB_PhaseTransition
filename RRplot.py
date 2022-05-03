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


data_constp = np.genfromtxt('Cornell22_const_P.csv', skip_header=1, delimiter=',')

p_constp = data_constp[:,0]
T_IC_constp = data_constp[:,1]
T_HEC_constp = data_constp[:,2]

t_HEC_constp = T_HEC_constp/h.Tc_mK(p_constp)
t_IC_constp = T_IC_constp/h.Tc_mK(p_constp)

Rcrit_HEC = h.critical_radius(t_HEC_constp, p_constp)
Rcrit_IC = h.critical_radius(t_IC_constp, p_constp)

Rcrit_min =  []

N=51

for p in p_constp:
    t_path = np.linspace(h.t_AB(p), 1.0 - 1/N, N)
    Rcrit_path = h.critical_radius(t_path, p)
    Rcrit_min.append(np.min(Rcrit_path))

Rcrit_min = np.array(Rcrit_min)

plt.scatter(Rcrit_min, Rcrit_HEC, marker='<', s=20, c='xkcd:pink', alpha=0.5, edgecolors='k', linewidth=1, label="HEC")
plt.scatter(Rcrit_min, Rcrit_IC, marker='<', s=15, c='xkcd:lightblue', alpha=0.5, edgecolors='b', linewidth=1, label="IC")

plt.xlim(0,4e4)
plt.ylim(0,3e4)
plt.grid()
plt.xlabel("Minimum AB interface curvature radius, A phase [nm]")
plt.ylabel(r"AB interface curvature radius at transition [nm]")
plt.legend(loc='upper left')