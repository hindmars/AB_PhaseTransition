#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:23:10 2022

He3 phase boundary solutions

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

# p = 25.5
p = 22
t = h.t_AB(p)

def get_and_plot_wall(t,p):

    A, pot, gr = hw.krylov_bubble(t,p, gr_pars=(500,250), dim=1)
    # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
    
    eden, eden_grad, eden_pot = hw.energy_density(A, pot, gr)
    
    sigma_AB = hw.energy(A,pot,gr)[0]*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    
    # Plotting
    fig, ax = plt.subplots(2,1, sharex='col')
    # gridspec_kw={'hspace': 0, 'wspace': 0}
    x = (gr.x - max(gr.x)/2)* h.xi(0,p)/h.xi(t,p)
    xmin = -20
    xmax = 20
    ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm()))
    
    ax[0].set_ylabel(r'$e/f_B$')
    ax[0].grid()
    ax[0].set_xlim(xmin, xmax)
    ax[0].set_title(r'AB boundary, p={:.1f} bar: $\sigma_{{AB}}/\xi_{{\rm GL}}(T)|f_B(T)| = {:.2f}$'.format(p, sigma_AB))
    
    norm =  np.sqrt(3)/h.delta_B_norm(t, p) 
    
    ax[1].plot(x, norm*A[:, 0, 0].real, label=r'${\rm Re}(A_{xx})$')
    ax[1].plot(x, norm*A[:, 1, 1].real, label=r'${\rm Re}(A_{yy})$')
    ax[1].plot(x, norm*A[:, 2, 2].real, label=r'${\rm Re}(A_{zz})$')
    ax[1].plot(x, norm*A[:, 0, 1].imag, label=r'${\rm Im}(A_{xy})$')
    ax[1].plot(x, norm*A[:, 1, 0].imag, label=r'${\rm Im}(A_{yx})$')
    
    ax[1].set_xlim(xmin, xmax)
    ax[1].legend()
    ax[1].grid()
    ax[1].set_ylabel(r'$A \sqrt{3}/\Delta_B(T,p)$')
    ax[1].set_xlabel(r'$x/\xi_{\rm GL}(T)$')

    return (A,pot,gr), ax

Apg_kry, ax = get_and_plot_wall(t,p)

sigma_AB = hw.energy(*Apg_kry)[0]*h.xi(0,p)/(abs(Apg_kry[1].mat_pars.f_B_norm())*h.xi(t,p))

print('Surface energy:', sigma_AB)
