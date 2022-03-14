#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 13:54:21 2022

He3 AB boundary - reproducing Schopol 1987

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

al = -1.0

# Schopol values 1
beta_hat_arr_1 = np.array([-0.288, 0.513, 0.504, 0.464, -0.643])
# Schopol values 1
beta_hat_arr_2 = np.array([-0.230, 0.461, 0.434, 0.384, -0.511])

A, pot, gr = hw.relax([0,10], al, beta_hat_arr_1, gr_pars=(500,20), dim=1)
# A, pot, gr = hw.krylov_bubble(al, beta_hat_arr_1,  gr_pars=(400,100), dim=1)
eden, eden_grad, eden_pot = hw.energy_density(A, pot, gr)

xi_ratio = h.xiGL_const/pot.mat_pars.xi()/(np.sqrt(5/3)) # Osheroff Cross definition

sigma_AB = hw.energy(A, pot, gr)[0]/abs(pot.mat_pars.f_B_norm())*(xi_ratio)

print('Surface energy:', sigma_AB)

x = (gr.x - max(gr.x)/2)* xi_ratio


fig, ax = plt.subplots(2,1, sharex='col')

# Plotting
xmin = -10
xmax = 10
ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm())/np.sqrt(5/3))

ax[0].set_ylabel(r'$e/f_B\sqrt{5/3}$')
ax[0].grid()
ax[0].set_xlim(xmin, xmax)
ax[0].set_title(r'AB boundary, $\sigma_{{AB}}/\xi_{{\rm GL}}^{{\rm OC}}(T)|f_B(T)| = {:.2f}$, relax method'.format(sigma_AB))

norm =  np.sqrt(3)/pot.mat_pars.delta_B_norm() 

ax[1].plot(x, norm*A[:, 0, 0].real, label=r'${\rm Re}(A_{xx})$')
ax[1].plot(x, norm*A[:, 1, 1].real, label=r'${\rm Re}(A_{yy})$')
ax[1].plot(x, norm*A[:, 2, 2].real, label=r'${\rm Re}(A_{zz})$')
ax[1].plot(x, norm*A[:, 0, 1].imag, label=r'${\rm Im}(A_{xy})$')
ax[1].plot(x, norm*A[:, 1, 0].imag, label=r'${\rm Im}(A_{yx})$')

ax[1].set_xlim(xmin, xmax)
ax[1].legend()
ax[1].grid()
ax[1].set_ylabel(r'$A \sqrt{3}/\Delta_B(T,p)$')
ax[1].set_xlabel(r'$x/\xi_{\rm GL}^{\rm OC}(T)$')
