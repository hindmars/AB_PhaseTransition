#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 14:41:17 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
from scipy.optimize import minimize
import scipy.optimize.nonlin
from scipy.integrate import solve_ivp

def F(tau, A_vec, alpha, beta_arr):
    """
    "Force" for gradient flow in He3 potential

    Parameters
    ----------
    tau : float
        ODE independent variable.
    A_vec : Complex array of shape (9,)
        Vectorised order parameter.

    Returns
    -------
    Complex array of shape (9,)
        Vectorised force.

    """
    A_mat = A_vec.reshape(3, 3)
    force = - h.dU_dA(A_mat, alpha, beta_arr)
    return force.reshape((9,))


def grad_flow(A_init, t, p, tau_max=25, n_tau=100):

    alpha = h.alpha_norm(t)
    beta_arr = h.beta_norm_asarray(t, p)
    
    tau = np.linspace(0, tau_max, n_tau)
    sol = solve_ivp(F, [min(tau),max(tau)], A_init.reshape((9,)), t_eval=tau, args=(alpha, beta_arr))
    
    A_vec_arr = np.transpose(sol.y) 
    A = A_vec_arr.reshape(A_vec_arr.shape[0], 3, 3)
    
    return A, tau


p = 28
t = 0.7#/h.Tc_mK(p)

v, A_AB, U_AB = h.line_section("A", "B", t, p)
v, A_Aplanar, U_Aplanar = h.line_section("A", "planar", t, p)
v, A_AX, U_AX = h.line_section("A", h.D_low, t, p)

ind_max = len(v)
n_flows = 5
ind_arr = (ind_max * np.linspace(1/n_flows, 1 - 1/n_flows, n_flows)).astype('int')

for ind in ind_arr:

    A_init = A_AX[ind]
    v_init = v[ind]
    
    A, tau = grad_flow(A_init, t, p, tau_max=50)
    
    # print("A_init", A_init)
    # print("A", A[-1])
    
    UA = h.U(A, h.alpha_norm(t), h.beta_norm_asarray(t, p))
    
    units = 1
    
    
    
    print("Initial energy", UA[0])
    print("Final Energy:     ", UA[-1].real )
    
    d_e_list = [UA[-1].real - h.f_phase_norm(t, p, phase) for phase in h.inert_phases]
    
    d_e_dict = dict(zip(h.inert_phases, d_e_list))
    
    d_e_dict_sorted = sorted(d_e_dict.items(), key = lambda kv : abs(kv[1]))
    
    print("nearest inert phase, energy diff", d_e_dict_sorted[0])
    # plt.plot(tau[0:UA.shape[0]], UA, label=r'$\phi_{{\rm i}}/\Delta_{{\rm wc}} = {:.2}$ ($\to$ {})'.format(v_init, d_e_dict_sorted[0][0]))
    plt.plot(h.norm(A) - h.delta_B_norm(t,p), UA, label=r'$\phi_{{\rm i}}/\Delta_{{\rm wc}} = {:.2}$ ($\to$ {})'.format(v_init, d_e_dict_sorted[0][0]))

# plt.xlabel(r'$\tau$')
plt.xlabel(r'$[\sqrt{{\rm tr} A A^\dagger} - \Delta_B]/(k_B T_c)$')
plt.ylabel(r'$f / f_0 $ ($T/T_c = 0.8$, $p = 26$ bar)')
plt.xlim(-0.5,0.5)
plt.ylim(-2,-1)
plt.plot([h.delta_A_norm(t,p) - h.delta_B_norm(t,p)]*2, [-2,-1],'k--', label="A phase gap")
plt.plot([-0.5, 0.5], [h.f_B_norm(t, p)]*2, 'k', label=r'$f_B/f_0$')

plt.title('Flows from minimum barrier path away from A phase')
plt.legend()
plt.grid()

