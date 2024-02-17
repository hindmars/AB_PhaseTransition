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


def bubble_action(*Apg):
    # A = Apg[0]
    pot = Apg[1]
    t = pot.mat_pars.t
    p = pot.mat_pars.p
    
    gr = Apg[2]
    x = np.ones_like(gr.x)
    
    A_phase = np.multiply.outer(x, h.D_A*h.delta_A_norm(t,p))
    
    E_bub = hw.energy(*Apg)[0]
    E_A = hw.energy(A_phase, pot, gr)[0]
    
    E_AB = E_bub - E_A
    
    S_AB =  E_AB * h.f_scale(p)*h.xi(0,p)**3/(h.kB * t * h.Tc_mK(p)*1e-3)

    return S_AB


def plot_bubble(A, pot, gr, plot_gap=False, loc='best'):
    
    eden, eden_grad, eden_pot = hw.energy_density(A, pot, gr)
    S_AB = bubble_action(A, pot, gr)
    xiGL = h.xi(pot.mat_pars.t, pot.mat_pars.p)

    # Plotting
    fig, ax = plt.subplots(2,1, sharex='col')
    # gridspec_kw={'hspace': 0, 'wspace': 0}
    # x = (gr.x - max(gr.x)/2)* h.xi(0,p)/h.xi(t,p)
    x = gr.x * h.xi(0, pot.mat_pars.p)
    xmin = min(x)
    xmax = max(x)
    ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm()) + 1)
    
    ax[0].set_ylabel(r'$e/f_B$ + 1')
    ax[0].grid()
    ax[0].set_xlim(xmin, xmax)
    title_string = r'AB bubble, p={:.1f} bar: $T/T_{{\rm c}} = {:.2f}$, $\xi_{{\rm GL}}(T) = {:.1f}$nm'.format(p, t, xiGL)
    title_string += '\n' + r' Barrier exponent: $E_{{\rm cb}}/T = {:.0f}$'.format(S_AB)
    ax[0].set_title(title_string)
    
    norm =  np.sqrt(3)/h.delta_B_norm(t, p) 
    
    ax[1].plot(x, norm*A[:, 0, 0].real, label=r'${\rm Re}(A_{xx})$')
    ax[1].plot(x, norm*A[:, 1, 1].real, label=r'${\rm Re}(A_{yy})$')
    ax[1].plot(x, norm*A[:, 2, 2].real, label=r'${\rm Re}(A_{zz})$')
    ax[1].plot(x, norm*A[:, 0, 1].imag, label=r'${\rm Im}(A_{xy})$')
    ax[1].plot(x, norm*A[:, 1, 0].imag, label=r'${\rm Im}(A_{yx})$')
    
    if plot_gap:
        ax[1].plot(x, norm * h.norm(A), 'k--', label=r'$||A||$')
    
    ax[1].set_xlim(xmin, xmax)
    ax[1].legend(loc=loc)
    ax[1].grid()
    ax[1].set_ylabel(r'$A \sqrt{3}/\Delta_B(T,p)$')
    ax[1].set_xlabel(r'$r$/nm')
    
    return ax


def get_and_plot_bubble(t, p, gr_pars=(500, 50), plot_gap=False, loc='best', maxiter=200):

    A, pot, gr = hw.krylov_bubble(t, p, gr_pars=gr_pars, dim=3, maxiter=maxiter)
    # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
        
    # sigma_AB = hw.energy(A,pot,gr)[0]*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    
    ax = plot_bubble(A, pot, gr, plot_gap, loc)


    return (A,pot,gr), ax

#%%

# h.set_default('DEFAULT_T_SCALE', 'Greywall')

p = 5.5
# p = 32
# t = 0.7*h.t_AB(p)
t = 0.98

Apg_kry, ax = get_and_plot_bubble(t, p, (500,100), plot_gap=True, maxiter=400)
S_AB = bubble_action(*Apg_kry)

print(f'Bubble energy/kBT (t = {t:.2f}, p = {p:.2f}): {S_AB:.0f}')

#%%    

plot_bubble(*Apg_kry, True, loc='lower left')