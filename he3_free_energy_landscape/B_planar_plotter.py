#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:25:44 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


p = 32
t = 0.7#/h.Tc_mK(p)


def line_section(X, D, t, p, n=500):
    v = np.linspace(0,1,n)*2
  
    if isinstance(X, str):
        X = h.D_dict[X]
        if isinstance(D, str):
            D = h.D_dict[D] - X

    A_XD = np.multiply.outer( np.ones_like(v) , X) + np.multiply.outer( v , D)
    U_XD = h.U(A_XD, (t-1), h.beta_norm_asarray(t, p) )
    
    return v, A_XD, U_XD

v = np.linspace(0,1,1000)*2

delta_A = h.delta_A_norm(t,p)
delta_B = h.delta_B_norm(t,p)
delta_planar = h.delta_planar_norm(t,p)
delta_polar = h.delta_polar_norm(t,p)





# The following is likely to be "closer" to the B phase than the canonical D_A
D_AX = np.matmul(h.O_xz, h.D_A)
D0 = (np.outer(h.e[+1], h.e[0]) + np.outer(h.e[-1], h.e[0]).conj())*0.5

A_A = np.multiply.outer( v , D_AX * delta_A)
A_B = np.multiply.outer( v , h.D_B * delta_B)


A_AB = np.multiply.outer( (1 - v) , D_AX * delta_A) + np.multiply.outer( v , h.D_B * delta_B)

A_A0 = np.multiply.outer( np.ones_like(v) , h.D_dict["A"] * delta_A) + np.multiply.outer( v , D0 * delta_A)

A_Agamma = np.multiply.outer( (1 - v) , D_AX * delta_A) + np.multiply.outer( v , h.D_dict["gamma"] * delta_B)
A_Apolar = np.multiply.outer( (1 - v) , D_AX * delta_A) + np.multiply.outer( v , h.D_dict["polar"] * delta_polar)

A_Aplanar = np.multiply.outer( (1 - v) , D_AX * delta_A) + np.multiply.outer( v , h.D_planar * delta_planar)
A_Bplanar = np.multiply.outer( (1 - v) , h.D_B * delta_B) + np.multiply.outer( v , h.D_planar * delta_planar)


U_Agamma = h.U(A_Agamma, (t-1), h.beta_norm_asarray(t, p) )
U_Apolar = h.U(A_Apolar, (t-1), h.beta_norm_asarray(t, p) )

U_A = h.U(A_A, (t-1), h.beta_norm_asarray(t, p) )
U_B = h.U(A_B, (t-1), h.beta_norm_asarray(t, p) )
U_AB = h.U(A_AB, (t-1), h.beta_norm_asarray(t, p) )
U_A0 = h.U(A_A0, (t-1), h.beta_norm_asarray(t, p) )
U_Aplanar = h.U(A_Aplanar, (t-1), h.beta_norm_asarray(t, p) )
U_Bplanar = h.U(A_Bplanar, (t-1), h.beta_norm_asarray(t, p) )

# v, A, U = line_section("A", "B", t, p)

# plt.plot(v, U)

#Jm3 = h.f_scale(p)*1e27
plt.plot(v,  U_AB , label=r'$X = A, D= A \to B$')
plt.plot(v,  U_Aplanar , label=r'$X = A, D= A \to$ planar')
# plt.plot(v,  U_Bplanar , label='X = planar, Y=B')
plt.plot(v,  U_A0 , label=r'$X = A, D= (+,0)$')

# plt.plot(v,  U_AB , label='X = B, Y=A')
# plt.plot(v,  U_Aplanar , label='X = planar, Y=A')
# plt.plot(v,  U_Bplanar , label='X = planar, Y=B')
# plt.plot(v,  U_Agamma, label=r'X = $\gamma$ (?)')
# plt.plot(v,  U_Apolar, label=r'X = polar')
# plt.plot(v,  U_AX , label=r'X = ??')

plt.xlabel(r'$\phi/\phi_b$')
plt.ylabel(r'$f / [\frac{{1}}{{3}} N(0) (k_BT_c)^{{2}} ]$ ($T/T_c = {}$, $p = {}$ bar)'.format(t,p))
# plt.ylabel(r'$f$ / ( J/m3 ), ($T/T_c = {:.2f}$, $p = {:}$ bar)'.format(t,p))
# plt.xlim(0,max(v))
# plt.xlim(0,0.1)
plt.ylim(-2,0)
# plt.ylim(-15,10)
plt.legend()
plt.grid()

plt.savefig("U_slices_from_A.pdf")