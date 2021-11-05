#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:43:15 2021

@author: hindmars
"""

"""
Minimal script for demonstrating use of critical_bubble_minimal module to 
find a critical bubble.

Main script is critical_bubble_minimal.krylov_bubble

It returns the solution by either 
- specifying the order parameter value at the maximum 
(for a quartic potential with quartic coupling lambda = 1.0 and 
broken phase order parameter phi_b = 1).
- specifying a quartic_paramset object, and a temperature at which to evaluate it.

"""


import critical_bubble_minimal as cbm
import matplotlib.pyplot as plt


# First, find the bubble solutions the simple way and find its action:
phi_m = 0.4
phi_1, pot_1, gr_1 = cbm.krylov_bubble(phi_m, display=True)

print("Actions components (total, gradient, potential) for bubble 1 with phi_m = {}".format(phi_m))
print(cbm.action(phi_1, pot_1, gr_1))




# Now, via a param_set, chosen "at random"
T0 = 0.5
D = 1
E = 2
lam = 1
T = 0.6
qpp = cbm.quartic_paramset(T0, D, E, lam)
phi_2, pot_2, gr_2 = cbm.krylov_bubble(qpp, T=T, display=True)

print("Action components (total, gradient, potential) for bubble 2 with:")
print("T0={}\nD={}\nE={}\nlambda={}\nT={}".format(T0,D,E,lam,T))
print(cbm.action(phi_2, pot_2, gr_2))

