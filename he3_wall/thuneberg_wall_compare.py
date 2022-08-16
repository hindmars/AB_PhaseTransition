#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 11:43:44 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

# h.set_default("DEFAULT_SC_CORRS", "WS15")


p_array = np.linspace(22,34,13)

surface_energy_us = []
surface_energy_thune = []

N, L_xi = (500, 10)

for p in p_array:
    t = h.tAB(p)
    A, pot, gr = hw.krylov_bubble(t,p, gr_pars=(N,L_xi*h.xi(t,p)), dim=1)
    # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
    
    eden, eden_grad, eden_pot = hw.energy_density(A, pot, gr)
    
    sigma_AB = hw.surface_energy(A,pot,gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    print(p, sigma_AB)
    surface_energy_us.append(sigma_AB)
    surface_energy_thune.append(hw.thuneberg_formula(t, p)/(-h.f_B_norm(t, p)*h.xi(t,p)))
    
#%%

surface_energy_us = np.array(surface_energy_us)
surface_energy_thune = np.array(surface_energy_thune)

plt.figure(figsize=(5,3))

plt.plot(p_array, surface_energy_us, label=r'Newton-Krylov optimisation, $N={}$, $L/\xi={}$'.format(N,L_xi))
plt.plot(p_array, surface_energy_thune, label='Thuneberg 1991')
plt.xlabel(r'$p$/bar')
plt.ylabel(r'$\sigma_{AB}/|f_B|\xi_{GL}(T)$')
plt.grid(True)
plt.xlim(22,34)
plt.ylim(0.9,1.0)
plt.title("AB interface surface surface energy")
plt.legend()
