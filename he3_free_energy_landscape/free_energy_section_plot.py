#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:25:44 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


def plot_free_energy_sections(t, p, savefig=False):

    fig, ax = plt.subplots()
    v, A_AB, U_AB = h.line_section("A", "B", t, p)
    v, A_AX, U_AX = h.line_section("A", h.D_low, t, p)
    v, A_Bplanar, U_Bplanar = h.line_section("B", "planar", t, p)
    
    # Construct the Osheroff-Cross 1977 path from A phase to planar phase
    
    A_Axz = np.matmul(h.D_dict['A'],h.O_yz)*h.delta_A_norm(t, p)
    A_planar = h.D_dict["planar"]*h.delta_planar_norm(t, p)
    v, A_Aplanar, U_Aplanar = h.line_section(A_Axz, A_planar, t, p, path='OC77')
    
    
    # Normalise free energies
    
    fA = h.f_phase_norm(t, p, "A")
    diff_fAB = np.abs(fA - h.f_phase_norm(t,p,"B"))
    scale = 1/diff_fAB
    
    
    def phi(X, X0=None):
        if X0 is None:
            X0 = X[0]
        return h.distance(X, X0)/h.delta_B_norm(t, p)
    
    
    ax.plot(phi(A_AB),  scale*(U_AB - fA), label=r'$A \to$ $B$ (straight path)')
    ax.plot(phi(A_Aplanar),  scale*(U_Aplanar - fA), label=r'$A \to$ planar (OC77 path)')
    ax.plot(phi(A_Bplanar),  scale*(U_Bplanar - fA), label=r'B $\to$ planar')
    ax.plot(phi(A_AX),  scale*(U_AX - fA), label=r'$A$ in direction $D_{\rm low}$')
    
    ind1 = np.argmin(h.distance(A_Aplanar, A_planar))
    ax.scatter(phi(A_Aplanar)[ind1], scale*(U_Aplanar - fA)[ind1], marker='o', color='k', label='planar phase' )
    
    ind2 = np.argmin(h.distance(A_Bplanar, A_planar))
    ax.scatter(phi(A_Bplanar)[ind2], scale*(U_Bplanar - fA)[ind2], marker='o', color='k')
    
    ax.set_xlabel(r'$||A - A_{0}||/\Delta_{\rm B}$')
    ax.set_ylabel(r'$(f - f_A) / |f_A - f_B|$')
    
    # Direct AB barrier is likly to be the highest and set y axis scale
    indB = np.argmin(np.abs(v - 1 ))
    fAB_barrier = np.max(U_AB[0:indB] - fA)
    # print(fAB_barrier)

    ax.set_xlim(0,max(v))
    ax.set_ylim(-1, np.ceil(scale*fAB_barrier / 2) * 2)
    
    ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Free energy on paths in order parameter space\n' + r'($T = {:.2f}$ mK, $p = {:.1f}$ bar)'.format(t*h.Tc_mK(p), p))
    
    if savefig:
        plt.savefig("free_energy_sections_t{:.2f}_p{:.1f}.pdf".format(t, p))

    return fig, ax