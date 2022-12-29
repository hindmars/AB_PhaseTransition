#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 12:04:16 2022

Plots Tc and TAB in Greywall and PLTS scales.

Also checks method for changing default temperature scale applied in module.

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

p = np.linspace(0, 34, 100)


# h.DEFAULT_T_SCALE="Greywall"
h.set_default('DEFAULT_T_SCALE', "Greywall")
Tc_mK_G = h.Tc_mK_expt(p)
TAB_mK_G = h.TAB_mK_expt(p)

# h.DEFAULT_T_SCALE="PLTS"
h.set_default('DEFAULT_T_SCALE', "PLTS")
Tc_mK_PLTS = h.Tc_mK_expt(p)
TAB_mK_PLTS = h.TAB_mK_expt(p)
# h.DEFAULT_T_SCALE="Greywall"

line_G = plt.plot(Tc_mK_G, p, label="$T_c$ Greywall")
line_PLTS = plt.plot(Tc_mK_PLTS, p, label="$T_c$ PLTS")

plt.plot(TAB_mK_G, p, label="$T_{AB}$ Greywall", ls='--', color=line_G[0]._color)
plt.plot(TAB_mK_PLTS, p, label="$T_{AB}$ PLTS", ls='--', color=line_PLTS[0]._color)


plt.grid()
plt.legend()
plt.xlabel(r'$T$/mK')
plt.ylabel(r'$p$/bar')