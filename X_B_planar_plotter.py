#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:25:44 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


#25,5	2,0333	2,04084

# p = 28.0
# t = 0.7 #0.1*h.Tc_mK(p)

p = 25.5
# t = 2.0333/h.Tc_mK(p)
t = 0#ß2.04084/h.Tc_mK(p)


v, A_AB, U_AB = h.line_section("A", "B", t, p)
v, A_Aplanar, U_Aplanar = h.line_section("A", "planar", t, p)
v, A_AX, U_AX = h.line_section("A", h.D_low, t, p)



# Jm3 = h.f_scale(p)*1e27
# Jm3 = -1/h.f_phase_norm(t,p,"B")
diff_fAB = h.f_phase_norm(t,p,"A") - h.f_phase_norm(t,p,"B")
scale = 1/diff_fAB


plt.plot(v,  scale*(U_AB - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = A \to B$')
plt.plot(v,  scale*(U_Aplanar - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = A \to$ planar')
plt.plot(v,  scale*(U_AX - h.f_phase_norm(t,p,"A")), label=r'$X = A, D = D_{\rm low}$')

ind_max = len(v)
n_flows = 5
ind_arr = (ind_max * np.linspace(1/n_flows, 1 - 1/n_flows, n_flows)).astype('int')


for ind in ind_arr:
    plt.scatter(v[ind], scale*(U_AX[ind] - h.f_phase_norm(t,p,"A")), marker='o' )

plt.xlabel(r'$\phi/\Delta_{\rm wc}$')
# plt.ylabel(r'$f / [\frac{1}{3} N(0) (k_BT_c)^{2} ]$ ($T/T_c = 0.8$, $p = 26$ bar)')
plt.ylabel(r'$(f - f_A) / (f_A - f_B)$ ($T = {:.2f}$ mK, $p = {:.1f}$ bar)'.format(t*h.Tc_mK(p), p))
# plt.ylabel(r'$(f - f_A)$ / ( J/m3 ), ($T/T_c = {:.2f}$, $p = {:}$ bar)'.format(t,p))
plt.xlim(0,max(v))
plt.ylim(-1,5)
# plt.ylim(-2,2)
plt.legend()
plt.grid()

plt.savefig("free_energy_sections.pdf")
