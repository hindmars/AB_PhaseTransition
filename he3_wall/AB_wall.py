#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:23:10 2022

He3 phase boundary solutions

@author: hindmars
"""

import he3_tools as h
import he3_wall as hw

# h.set_default("DEFAULT_SC_CORRS", "WS15")
# p = 25.5
p = 29
t = h.t_AB(p)

N, L_xi = (500, 25)

savefig=False

def get_and_plot_wall(t,p, w):

    (A, pot, gr), sigma_bw, sigma_tot = hw.get_wall(t, p, w)
    # A, pot, gr = hw.krylov_bubble(h.alpha_norm(t), h.beta_norm_asarray(t,p), gr_pars=(500,125), dim=1)
    
    # Plotting
    ax = hw.plot_wall(A, pot, gr)

    return (A,pot,gr), ax


Apg_kry, ax = get_and_plot_wall(t, p, L_xi * h.xi(t,p))

# sigma_AB = hw.energy(*Apg_kry)[0]*h.xi(0,p)/(abs(Apg_kry[1].mat_pars.f_B_norm())*h.xi(t,p))

if savefig:
    file_name = 'AB_wall_t={:.2f}_p={:.1f}.pdf'.format(t,p)
    print('print to figure',file_name)
    ax[0].get_figure().savefig(file_name)

