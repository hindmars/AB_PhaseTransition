#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 14:41:17 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
from scipy.optimize import newton_krylov, minimize
import scipy.optimize.nonlin


p = 28
t = 0.7#/h.Tc_mK(p)

# A_init = h.D_dict["planar"]*3
# A_init = h.D_dict["B"]*3

# A_init = 3*(np.random.rand(3,3)-0.5 + 1j*np.random.rand(3,3)-0.5j)

v, A_AB, U_AB = h.line_section("A", "B", t, p)
v, A_Aplanar, U_Aplanar = h.line_section("A", "planar", t, p)
v, A_AX, U_AX = h.line_section("A", h.D_low, t, p)

ind = 271

A_init = A_AX[ind]

# A_init = h.D_dict["A"] * h.delta_phase_norm(t, p, "A") + h.D_low * 

def F(A):
    return h.dU_dA(A, h.alpha_norm(t), h.beta_norm_asarray(t, p))


# def U_vec(v, alpha, beta_arr):
#     return h.U(arr2cmat(v), alpha, beta_arr)

# def cmat2arr(A):
#     v=np.zeros((18,))
#     v[0:9] = np.reshape(A.real, (9,))
#     v[9:18] = np.reshape(A.imag, (9,))
#     return v

# def arr2cmat(v):
#     u_re = np.reshape(v[0:9], (3,3))
#     u_im = np.reshape(v[9:18], (3,3))
#     A = u_re + 1j*u_im
#     return(A)    

try:
    A = newton_krylov(F, A_init)
except scipy.optimize.nonlin.NoConvergence as e:
    A = e.args[0]
    print('No Convergence')

# v_init = cmat2arr(A_init)

# res = minimize(U_vec, v_init, args=( h.alpha_norm(t), h.beta_norm_asarray(t, p)))

# A_min = arr2cmat(res.x)

print("A_init", A_init)
print("A", A)
# print("A_min", A_min)
# print(h.D_dict["planar"] * h.delta_phase_norm(t,p,"B"))
# print(h.D_dict["B"] * h.delta_phase_norm(t,p,"B"))

UA = h.U(A, h.alpha_norm(t), h.beta_norm_asarray(t, p))

Jm3 = 1

plt.plot(v,  Jm3*(U_AB - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = A \to B$')
plt.plot(v,  Jm3*(U_Aplanar - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = A \to$ planar')
plt.plot(v,  Jm3*(U_AX - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = D_{\rm low}$')

plt.scatter(v[ind], Jm3*(U_AX[ind] - h.f_phase_norm(t,p,"A")), marker='o' )

plt.xlabel(r'$\phi/\Delta_{\rm wc}$')
plt.ylabel(r'$f / [\frac{1}{3} N(0) (k_BT_c)^{2} ]$ ($T/T_c = 0.8$, $p = 26$ bar)')
# plt.ylabel(r'$(f - f_A)$ / ( J/m3 ), ($T/T_c = {:.2f}$, $p = {:}$ bar)'.format(t,p))
plt.xlim(0,max(v))
plt.ylim(-0.1,0.5)
# plt.ylim(-2,2)
plt.legend()
plt.grid()

print("Energy:     ", UA.real )
print("Energy diff:", UA.real - h.f_phase_norm(t,p,"A"))

for phase in h.inert_phases:
    print(phase, h.f_phase_norm(t,p,phase))
    
