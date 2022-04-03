#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 19:02:45 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

p = np.linspace(0,34,100)

# TAB_PLTS = h.TAB_poly_PLTS(p)
TAB_Gw = h.TAB_poly_Greywall(p - h.p_pcp_bar)
Tc_Gw = h.Tc_poly_Greywall(p)

# TAB_PLTS[p<h.p_pcp_bar] = np.nan
TAB_Gw[p<h.p_pcp_bar] = np.nan

plt.plot(Tc_Gw, p, 'k', label=r"$T_{c}$, Greywall 1986")
plt.plot(TAB_Gw, p, 'k--', label=r"$T_{AB}$, Greywall 1986")
# plt.plot(TAB_PLTS, p, label=r"$T_{AB}$, PLTS")


TAB_RWS = h.t_AB(p)*h.Tc_mK(p)
plt.plot(TAB_RWS, p, label=r"$T_{AB}$, RWS 2019")

h.DEFAULT_SC_CORRS = "Wiman_thesis"
TAB_Wth = h.t_AB(p)*h.Tc_mK(p)
plt.plot(TAB_Wth, p, label=r"$T_{AB}$, Wiman thesis")
h.DEFAULT_SC_CORRS = "RWS19"

plt.grid()
plt.legend()
plt.xlim(1.8, 2.5)
plt.ylim(21, 30)
# plt.ylim(0,34)
plt.xlabel(r"$T$/mK")
plt.ylabel(r"$p$/bar")


