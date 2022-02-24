#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 08:59:01 2022

Script for mapping out bubble action over a range of lambda/bar


@author: hindmars
"""

import critical_bubble_minimal as cbm
import numpy as np
import matplotlib.pyplot as plt


def get_psi_m(lam_bar):
    return (3 - np.sqrt(9-8*lam_bar))/(2*lam_bar)

# First, find the bubble solutions the simple way and find its action:
lam_bar_arr = np.linspace(0.01, 0.99, 49)
action_scaled_arr = np.zeros_like(lam_bar_arr)

phi_list = []
pot_list = []

for n, lam_bar in enumerate(lam_bar_arr):
    phi, pot, gr = cbm.krylov_bubble(lam_bar, gr_pars=(2000,200), dim=3)
    scale = pot.mu/pot.lam**1.5
    action_scaled_arr[n] = cbm.action(phi, pot, gr)[0]/scale
    phi_list.append(phi)
    pot_list.append(pot)

    print("Action lambda_bar = {}:  {}".format(lam_bar, action_scaled_arr[n]))

fig, ax = plt.subplots()
ax.semilogy(lam_bar_arr, action_scaled_arr, '.')
ax.set_xlabel(r'$\bar\lambda$')
ax.set_ylabel(r'$E_c \lambda^{3/2}/ \delta$')
ax.set_xlim(0,1)
ax.grid()

fig.savefig("bubble_energy.pdf")