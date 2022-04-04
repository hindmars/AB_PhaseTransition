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

p = np.linspace(0,34,500)
plt.figure()

Tc_exp = h.Tc_mK_expt(p)
TAB_exp = h.TAB_mK_expt(p)

plt.plot(Tc_exp, p, 'k', label=r"$T_{c}$, Greywall 1986")
plt.plot(TAB_exp, p, 'k--', label=r"$T_{AB}$, Greywall 1986")

TAB_RWS = h.t_AB(p)*h.Tc_mK(p)
plt.plot(TAB_RWS, p, label=r"$T_{AB}$, RWS 2019")


h.DEFAULT_SC_CORRS = "Wiman-thesis"

TAB_Wth = h.t_AB(p)*h.Tc_mK(p)
plt.plot(TAB_Wth, p, label=r"$T_{AB}$, Wiman thesis")

h.DEFAULT_SC_CORRS = "RWS19"
p_coarse = p[::5]
# TAB_RWS_adj = h.t_AB(p_coarse, sc_adjust=True)*h.Tc_mK(p_coarse)
h.DEFAULT_SC_ADJUST=True
TAB_RWS_adj = h.t_AB(p_coarse)*h.Tc_mK(p_coarse)
h.DEFAULT_SC_ADJUST=False
plt.plot(TAB_RWS_adj, p_coarse, 'x', label=r"$T_{AB}$, RWS 2019, adjusted")



plt.grid()
plt.legend(loc='upper right')
plt.xlim(1.8, 2.5)
plt.ylim(21, 34)
# plt.ylim(0,34)
plt.xlabel(r"$T$/mK")
plt.ylabel(r"$p$/bar")
plt.title("Temperature scale " + h.DEFAULT_T_SCALE)


