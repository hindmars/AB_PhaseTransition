#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 15:54:05 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import he3_wall as hw

def plot_confine(A, pot, gr):
    
    t = pot.mat_pars.t
    p = pot.mat_pars.p
    eden, eden_grad, eden_pot = hw.energy_density(A, pot, gr)
    sigma_tot = np.trapz(eden, gr.x)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    sigma_bw = hw.surface_energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))

    fig, ax = plt.subplots(2,1, sharex='col')
    # gridspec_kw={'hspace': 0, 'wspace': 0}
    x =  gr.x * h.xi(0,p)/h.xi(t,p)
    xmin = min(x)
    xmax = max(x)
    
    w_mu = gr.R * h.xi(0,p)/1000
    
    ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm()) + 1)
    
    ax[0].set_ylabel(r'$e/|f_B| + 1$')
    ax[0].grid()
    ax[0].set_xlim(xmin, xmax)
    tstring1 = r'p={:.1f} bar, $T={:.2f}$ mK, $w={:.2f} \;\mu$m'.format(p, t*h.Tc_mK(p), w_mu)
    tstring2 = r'$\sigma_{{bw}}/\xi_{{\rm GL}}(T)|f_B(T)| = {:.3f}$, '.format(sigma_bw)
    tstring3 = r'$\sigma_{{tot}}/\xi_{{\rm GL}}(T)|f_B(T)| = {:.3f}$'.format(sigma_tot)
    ax[0].set_title(tstring1 + '\n' + tstring2 + tstring3 )
    
    norm =  np.sqrt(3)/h.delta_B_norm(t, p) 
    
    ax[1].plot(x, norm*A[:, 0, 0].real, label=r'${\rm Re}(A_{xx})$')
    ax[1].plot(x, norm*A[:, 1, 1].real, label=r'${\rm Re}(A_{yy})$')
    ax[1].plot(x, norm*A[:, 2, 2].real, label=r'${\rm Re}(A_{zz})$')
    ax[1].plot(x, norm*A[:, 1, 2].imag, label=r'${\rm Im}(A_{yz})$')
    # ax[1].plot(x, norm*A[:, 2, 1].imag, label=r'${\rm Im}(A_{zy})$')
    ax[1].plot(x, norm*A[:, 1, 0].real, label=r'${\rm Re}(A_{yx})$')
    
    ax[1].set_xlim(xmin, xmax)
    ax[1].legend(loc='center')
    ax[1].grid()
    ax[1].set_ylabel(r'$A \sqrt{3}/\Delta_B(T,p)$')
    ax[1].set_xlabel(r'$x/\xi_{\rm GL}(T)$')

    return ax


def get_confine(t, p, w, wall_phase="Ay", bulk_phase="B",
                bcs = [hw.bc_min_pb, hw.bc_min_pb], T_list = None, 
                N=500, **kwargs):
    
    L = (w)/h.xi(0,p)
    
    A, pot, gr = hw.krylov_confine(t, p, wall_phase=wall_phase, bulk_phase=bulk_phase, 
                                   bcs=bcs, T_list=T_list, gr_pars=(N,L), dim=1, **kwargs)

    sigma_bw = hw.surface_energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    sigma_tot = hw.energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    # sigma_bw = sigma_bw/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    
    return (A, pot, gr), sigma_bw, sigma_tot


def get_and_plot_confine(t, p, w, wall_phase="Ay", bulk_phase="B",
                         bcs = [hw.bc_min_pb, hw.bc_min_pb], T_list = None, 
                         N=500, **kwargs):

    Apg, sigma_bw, sigma_tot = get_confine(t, p, w, wall_phase, bulk_phase, 
                                       bcs, T_list, 
                                       N, **kwargs)
    
    # Plotting
    ax = plot_confine(*Apg)

    return Apg, sigma_bw, sigma_tot, ax

# # p = 25.5
# p = 26
# t = 1.95/h.Tc_mK(p)
# w = 1100

# Apg_conf_BB, sigma_bw_BB, sigma_tot_BB, _ = get_and_plot_confine(t, p, w/2, "B", "B", bcs = [hw.bc_max_pb, hw.bc_neu])
# Apg_conf_AyB, sigma_bw_AyB, sigma_tot_AyB, _ = get_and_plot_confine(t, p, w/2, "Ay", "B", bcs = [hw.bc_max_pb, hw.bc_neu])
# Apg_conf_AyAy, sigma_bw_AyAy, sigma_tot_AyAy, _ = get_and_plot_confine(t, p, w/2, "Ay", "Ay", bcs = [hw.bc_max_pb, hw.bc_neu])

# print('AyAy: Pressure:', p, 'T:', t*h.Tc_mK(p), ', Surface energy:',  sigma_bw_AyAy, 'Total energy:',  sigma_tot_AyAy[0])
# print('AyB: Pressure:', p, 'T:', t*h.Tc_mK(p), ', Surface energy:',  sigma_bw_AyB, 'Total energy:',  sigma_tot_AyB[0])
# print('BB: Pressure:', p, 'T:', t*h.Tc_mK(p), ', Surface energy:',  sigma_bw_BB, 'Total energy:',  sigma_tot_BB[0])
