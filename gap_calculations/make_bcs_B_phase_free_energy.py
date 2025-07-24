#!/usr/bin/env, python3
#, -*-, coding:, utf-8, -*-
"""
Created on Wed Aug 10 08:59:45 2022

BCS gap calculator

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
import scipy.optimize as scopt
import scipy.special as scspe

import he3_tools as h

import bcs_tools as bcs

r""" The BCS equation for the free emergy $F$ is (Vollhardt and Woelfle sect 3.3.4)
$$
1 = - V_0 N(0) \int_0^{\epsilon_c} d\xi \frac{1}{E} 
    \tanh \left( \frac{E}{2k_B T} \right)
$$
where $\xi$ is the quasiparticle energy. $E = \sqrt{\xi^2 + \Delta^2}$.
$\epsilon_c = \hbar \omega_D$, and $\omega_D$ is the Debye frequency.

At the critical temperature $T_c$, the gap vanishes, and so
$$
1 = - V_0 N(0) \int_0^{\epsilon_c} d\xi \frac{1}{\xi} 
    \tanh \left( \frac{\xi}{2k_B T_c} \right).
$$
This gives the critical tempature in terms of the scattering potential, the 
density of states, and the cut-off (which can be viewed as part of the scattering 
potential).

In order to solve for the gap, we subtract and divide by $V_0N(0)$, to obtain
$$
\int_0^{\epsilon_c} d\xi 
\left[ \frac{1}{E} \tanh \left( \frac{E}{2k_B T} \right) - 
      \frac{1}{\xi} 
          \tanh \left( \frac{\xi}{2k_B T_c} \right)
      \right] = 0.
$$
The integral is convergent, and insensitive to the cut-off, which can be taken 
to infinity.  $x = \xi/k_BT_c$, $D = \Delta/k_BT_c$, and $t = T/T_c$, we have 
to solve the equation
$$
\Gamma(D,t) = \int_0^{\infty} d x \left[ \frac{1}{\sqrt{x^2 + D^2}} 
    \tanh \left( \frac{1}{2t}\sqrt{x^2 + D^2} \right) - 
      \frac{1}{x} 
          \tanh \left( \frac{x}{2} \right)
      \right] = 0 .
$$
This defines the function $D(t)$, and hence the BCS gap is
$$
\Delta_T = D(T/T_c) k_BT_c.
$$
A good fit to the function is 
$$
\Delta_T^\text{fit1} = \Delta_0 \left(1 - t^a\right)^{0.5}, \qquad a = 3.285.
$$
Not noticeably better is 
$$
\Delta_T^\text{fit2} = \Delta_0 \left(1 - t^a\right)^b, \qquad a = 3.485,\quad b = 0.541.
$$
Polynomials also work
$$
\Delta_T^\text{fit3} = \Delta_0 \left(a_0 + a_1 t + a_2 t^2 a_3 t^3)^{1/2}, 
\qquad a = (-2.06\times 10{-3},  3.14 -3.11,  0.961).
$$
"""




#%%

f_lim = 0.1
# t_arr = np.linspace(0.025,1.1,43)
t_arr = np.arange(0.025,1.1,0.025)
D_arr = np.linspace(0.0, 3.0, 31, endpoint=False)

#%% Get A phase gap

bcs_gap_B_phase_arr = np.array([float(scopt.fsolve(bcs.gap_eqn, 2, args=(t,))) for t in t_arr])


#%% Get B phase free energy

f_B_t_list = []

for t in t_arr:
    
    f_B_arr = np.zeros_like(D_arr)
    
    for n, D in enumerate(D_arr):
        f_B_arr[n] = bcs.free_energy(D, t)
    
    f_B_t_list.append(f_B_arr)

f_B_t_arr = np.array(f_B_t_list)


#%% plot gap

fig1, ax1 = plt.subplots()

ax1.plot(t_arr, bcs_gap_B_phase_arr**2)

ax1.set_xlabel(r'$T/T_c$')
ax1.set_ylabel(r'$(\Delta/k_BT_c)^2$')

ax1.set_xlim(0,1.1)
ax1.set_ylim(0,5)

ax1.grid(True)

ax1.set_title('BCS gap squared B phase')

#%% Plot free energy selected t

fig3, ax3 = plt.subplots()

for n in range(7,len(t_arr),8):
    
    f_B_arr = f_B_t_arr[n,:]
    
    popt, copt = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(-t, .1, .1))
    
    line = ax3.plot(D_arr, f_B_t_arr[n,:] - f_B_t_arr[n,0])

    ax3.plot(D_arr, bcs.poly_in_squares(D_arr, *popt), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.1f}$: '.format(t_arr[n]) + r'$a = ({:.2f}, {:.2f}, {:.3f})$'.format(*popt))


    print(popt[0]/3, popt[1]/9, popt[2]/27 )

ax3.set_xlabel(r'$\Delta/k_B T_c$')
ax3.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax3.set_xlim(0, 3)
ax3.set_ylim(-5, 5)

ax3.legend(fontsize='small')

ax3.set_title('BCS B-phase free energy (solid), fit to $a_0x^2 + a_1x^4 + a_2x^6$ (dash)' 
              + '\n'
              + r'Fitted $\beta_{{B}} = {:.5f}$  '.format(popt[1]/9)
              + r'$\beta_{{B}}^{{\rm BCS}} = {:.5f}$'.format(h.beta_B_norm(0, 0)), 
              fontsize='smaller')

ax3.grid(True)

#%% save gap. free energy

np.savez('bcs_gap_B_phase.npz', t=t_arr, bcs_gap_B_phase=bcs_gap_B_phase_arr)

np.savez('free_energy_B_phase.npz', t=t_arr, delta_B=D_arr, free_energy_B_phase=f_B_t_arr)

