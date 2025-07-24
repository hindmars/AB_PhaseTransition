#!/usr/bin/env, python3
#, -*-, coding:, utf-8, -*-
"""
Created on Wed Aug 10 08:59:45 2022

BCS gap calculator

@author: hindmars
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
import scipy.optimize as scopt
import scipy.special as scspe

from numpy.polynomial import Polynomial

import he3_tools as h
import he3_free_energy as h3fe
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

savefig = True
fig_dir = './figures/'

if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

#%% Load BCS calculations

bcs_gap_A_dict = np.load('bcs_gap_A_phase.npz')

bcs_f_A_dict = np.load('free_energy_A_phase.npz')

t_arr = bcs_f_A_dict['t']
D_arr = bcs_f_A_dict['delta_A']
f_A_t_arr = bcs_f_A_dict['free_energy_A_phase']

bcs_gap_A_phase_arr = bcs_gap_A_dict['bcs_gap_A_phase']

#%% Load RHUL gap data
# [1] Temperature, T [mK] 
# [2] Uncertainty in temperature, dT [mK]
# [3] Relative temperature, T/Tc0
# [4] Uncertainty in relative temperature
# [5] NMR precession frequency, f [Hz]
# [6] Frequency shift, f - fL [Hz]
# [7] Spatially averaged energy gap squared, <\Delta_A^2> / kB^2 = (IS_bulk/IS)*|f^2 - fL^2| [mK^2]
# [8] Uncertainty in squared energy gap [mK^2]

h.set_default('DEFAULT_T_SCALE', 'Greywall')
h.set_default('DEFAULT_SC_CORRS', 'Choi-interp')

data_dir = '/Users/hindmars/physics/He3/experiments/RHUL_gap/'

file_pattern = 'Run{:}_expr_gap_specular_{:}bar.dat'

run_list = [33]*3 + [34]*2
p_list = [0.17, 2.46, 5.50, 12.00, 21.00]
p_str_list = [f'{p:.2f}' for p in p_list]

color_dict = dict(zip(p_str_list, bcs.color_list))

data_dict = {}

for p_str, run in zip(p_str_list, run_list):
    data_dict[p_str] = np.loadtxt(data_dir + file_pattern.format(run, p_str), skiprows=20)

#%% Check Tc suppression with Tc from experiment

Tc_mK_RHUL_dict = {}
beta_A_RHUL_dict = {}

fig1, ax1 = plt.subplots()


def linear(x, *a):
    return a[0] + a[1]*x

for p_str in p_str_list: 

    p = float(p_str)
    Tc_mK =  h.Tc_mK(p)
    
    T = data_dict[p_str] [:, 0]
    Delta2_expt = data_dict[p_str] [:, 6]
    d_Delta2_expt = data_dict[p_str] [:, 7]

    just_below_Tc = (T < Tc_mK - 0.025) & (Tc_mK -0.05 < T)
    near_Tc = (np.abs(T - Tc_mK) < 0.05)
    
    popt, copt = scopt.curve_fit(linear, T[just_below_Tc], Delta2_expt[just_below_Tc], 
                                 p0=(-Tc_mK/bcs.bcs_const, 1/bcs.bcs_const))

    Tc_mK_est = -popt[0]/popt[1]
    beta_est = -1/popt[1]*Tc_mK_est/4

    ax1.errorbar(T[near_Tc]- Tc_mK, Delta2_expt[near_Tc], d_Delta2_expt[near_Tc], 
                 ls='', marker='.', c=color_dict[p_str], 
                 label=p_str + r': $T_c^{{\rm est}} = {:.2f}$ mK, $T_c^{{\rm G}} = {:.2f}$ mK'.format(Tc_mK_est, Tc_mK) )
    ax1.plot((T[just_below_Tc] - Tc_mK), linear(T[just_below_Tc], *popt), 'k--')

    print(p_str, Tc_mK, Tc_mK_est, beta_est, h.beta_A_norm(1,p))
    
    Tc_mK_RHUL_dict[p_str] = Tc_mK_est
    beta_A_RHUL_dict[p_str] = beta_est

ax1.set_xlabel(r'$(T - T_c^{\rm G})/$mK')
ax1.set_ylabel(r'$\Delta_A^2/k_{\rm B}^2$ [mK$^2$]')

ax1.set_title(r'$\Delta_A^2$ near $T_c$ with linear fits')
ax1.grid(True)
ax1.legend(fontsize='smaller')
ax1.set_xlim(-0.06,0.02)

fig1.tight_layout()

if savefig:
    fig1.savefig(fig_dir + 'gap2_A_near_Tc_fits.png')
    fig1.savefig(fig_dir + 'gap2_A_near_Tc_fits.pdf')

#%% Plot effective quartic parameter 

fig2, ax2 = plt.subplots()

def beta_bcs_approx(t, *a): 
    return 6/(1 + a[0]*(t) + (2 - a[0])*(t)**2)

alpha4 = t_arr - 1

beta4_gap = (1-t_arr)/bcs_gap_A_phase_arr**2/4

in_range = (t_arr < 1) & (t_arr > 0.4) 

popt_gap_beta_A, copt_gap_beta_A = scopt.curve_fit(beta_bcs_approx, 
                                           t_arr[in_range], beta4_gap[in_range]/h.beta_const, 
                                           p0=(1,))

line_beta_gap = ax2.plot(t_arr[t_arr<1], beta4_gap[t_arr<1]/h.beta_const, 
         # ls=bcs.ls_dict['long dash'], 
         c='r', 
         # label=r'$(1-t)/4\Delta^2_{\rm BCS}$')
         label=r'$\beta_A^{\rm BCS}$')

ax2.plot(t_arr, beta_bcs_approx(t_arr, *popt_gap_beta_A), 
         ls='--', c=line_beta_gap[0].get_c(), 
         label=r'$6/(1 + a_0t + (2-a_0)t^2)$' + '\n' + r'$a_0 = {:.2f}$'.format(*popt_gap_beta_A))

def beta_A_bcs_approx(t):
    return beta_bcs_approx(t, *popt_gap_beta_A)*h.beta_const


for p_str in p_str_list:
    p = float(p_str)
    T = data_dict[p_str] [:, 0]
    Delta2_expt = data_dict[p_str] [:, 6]
    d_Delta2_expt = data_dict[p_str] [:, 7]

    Tc_mK = Tc_mK_RHUL_dict[p_str]
    
    t_expt = T/Tc_mK
    beta4_expt = (1-t_expt)/(Delta2_expt/Tc_mK**2)/4
    d_beta4_expt = (1-t_expt)/((Delta2_expt - d_Delta2_expt)/Tc_mK**2)/4 - beta4_expt
    
    
    line_beta_expt = ax2.errorbar(t_expt[t_expt<1.0], beta4_expt[t_expt<1.0]/h.beta_const,
                                 d_beta4_expt[t_expt<1.0]/h.beta_const, 
                                 marker='.',
                                 markersize=1,
                                 ls='',
                                 c=color_dict[p_str],
             label=r'$p$/bar = ' + p_str)
    ax2.scatter(h.Tc_mK(p)/Tc_mK, h.beta_A_norm(1,p)/h.beta_const, 
                c=line_beta_expt[0].get_c())


ax2.set_title(r'Effective quartic parameter $\beta_A^{\rm eff} = (1 - t)/2 \Delta_A^2(t)$'
              + '\n'
              + 'Expt (errorbar), RWS19($T_c$) (dot)')

ax2.set_xlabel(r'$T/T_c$')
ax2.set_ylabel(r'$\beta_A^{\rm eff}(t)/\beta^{\rm wc}$')

ax2.set_xlim(0,1.1)
ax2.set_ylim(1.0, 4.0)

ax2.legend(fontsize='small')
ax2.grid(True)

fig2.tight_layout()

if savefig:
    fig2.savefig(fig_dir + 'beta_A_eff_t.png')
    fig2.savefig(fig_dir + 'beta_A_eff_t.pdf')


#%% Plot effective quartic parameter difference from BCS value

fig2a, ax2a = plt.subplots()

def beta_corr_fitfun(T, *a):
    val = np.zeros_like(T)
    for n, par in enumerate(a):
        val += par*T**(n)
    return val

def hei24_fun(t, p):
    beta_A_sc = h.beta_A_norm(1, p)
    beta_A_wc = h.beta_A_norm(0, p)
    
    d_beta = (beta_A_sc - beta_A_wc)/beta_A_sc
    
    return 1 - 0.37*d_beta*(1 - t**3)

alpha4 = t_arr - 1
beta4_gap = (1-t_arr)/bcs_gap_A_phase_arr**2/4

p0_dict = {}

p0_dict['0.17'] = (1.0,)
p0_dict['2.46'] = (1.0, 0.0)
p0_dict['5.50'] = (1.0, 0.0, 0.0)
p0_dict['12.00'] = (1.0, 0.0, 0.0)
p0_dict['21.00'] = (1.0, 0.0, 0.0)

popt_bc_dict = {}

derived_data = {}

for p_str in p_str_list:
    p = float(p_str)
    T = data_dict[p_str] [:, 0]
    Delta2_expt = data_dict[p_str] [:, 6]
    d_Delta2_expt = data_dict[p_str] [:, 7]

    Tc_mK = Tc_mK_RHUL_dict[p_str]
    # Tc_mK = h.Tc_mK(p)
    
    t_expt = T/Tc_mK
    beta4_expt = (1-t_expt)/(Delta2_expt/Tc_mK**2)/4
    d_beta4_expt = (1-t_expt)/((Delta2_expt - d_Delta2_expt)/Tc_mK**2)/4 - beta4_expt
        
    beta4_gap_interp = np.interp(t_expt, t_arr, beta4_gap)

    derived_data[p_str] = np.vstack((t_expt, beta4_gap_interp, beta4_expt, d_beta4_expt))

    beta4_diff = beta4_expt - beta4_gap_interp
    
    sc_corr_factor = h.beta_A_norm(1, p)/h.beta_A_norm(0, p)
    beta4_he124 = (sc_corr_factor*hei24_fun(t_arr, p))*beta4_gap
    beta4_hei24_diff = beta4_he124 - beta4_gap
    
    beta_A_wc = 2*h.beta_const
    beta4_diff_ratio = beta4_diff/beta_A_wc
    d_beta4_diff_ratio = np.abs(d_beta4_expt/beta_A_wc)
    beta4_hei24_diff_ratio = beta4_hei24_diff/beta_A_wc
    
    line_beta_expt = ax2a.errorbar(T[t_expt<1.0], 
                                  beta4_diff_ratio[t_expt<1.0], 
                                  d_beta4_diff_ratio[t_expt<1.0],
                              ls = '',
                              marker='.',
                              markersize=1,
                              c=color_dict[p_str],
                              label=p_str)

    d_beta_A_sc = h.beta_A_norm(1,p) - h.beta_A_norm(0,p)

    ax2a.plot(t_arr[t_arr<1.0]*Tc_mK, beta4_hei24_diff_ratio[t_arr<1.0], 'k--')

    ax2a.scatter(h.Tc_mK(p), d_beta_A_sc/beta_A_wc, 
                 c=color_dict[p_str])
    ax2a.scatter(Tc_mK_RHUL_dict[p_str], (beta_A_RHUL_dict[p_str] - h.beta_A_norm(0,p)) /beta_A_wc, 
                 c=line_beta_expt[0].get_c(),
                 marker='x')
    popt_bc, copt_bc = scopt.curve_fit(beta_corr_fitfun, T[t_expt<0.9], 
                                       beta4_diff_ratio[t_expt<0.9], 
                                       p0=p0_dict[p_str])
    
    print(p_str, popt_bc)

    popt_bc_dict[p_str] = popt_bc

    label_str = r'$\vec{a} = ('
    for n, a in enumerate(popt_bc):
        label_str += '{:.2f}, '.format(a)
        
    label_str += ')$'
    
    ax2a.plot(T[t_expt<1], beta_corr_fitfun(T[t_expt<1], *popt_bc), 
             c=color_dict[p_str],
             ls='--',
             label=label_str)


sc_params = h.get_setting('DEFAULT_SC_CORRS')

ax2a.set_title(r'Effective quartic parameter $\beta_A = (1 - t)/2 \Delta_A^2(t)$ expt (errorbar), '
              + '\n'
              + 'Heikkinen+24 (black dash), quadratic fit (dash), {:}($T_c$) (dot), fitted($T_c$) (cross)'.format(sc_params),
              fontsize='small')

ax2a.set_xlabel(r'$T$/mK')

ax2a.set_ylabel(r'$[\beta_A^{\rm expt}(t) - \beta_A^{\rm BCS}(t)]/\beta_A^{\rm wc}$')
ax2a.set_ylim(-0.5, 0.2)
ax2a.set_xlim(0.5, 2.5)

ax2a.grid(True)
ax2a.legend(title=r'$y = \sum_{n=0} a_n x^n$                   $p$/bar   ', 
            ncols=2, 
            fontsize='small', 
            alignment='right')

fig2a.tight_layout()

if savefig:
    fig2a.savefig(fig_dir + 'beta_A_eff_diff_T.png')
    fig2a.savefig(fig_dir + 'beta_A_eff_diff_T.pdf')


#%% Plot quartic aproximation gap and compare to data


fig3, ax3 = plt.subplots()

tr_AdagA = 2

for p_str in p_str_list[::-1]: 
    
    p = float(p_str)
    # Tc_mK =  h.Tc_mK(p)

    Tc_mK =  Tc_mK_RHUL_dict[p_str]
    sc_corr_factor = h.beta_A_norm(1, p)/h.beta_A_norm(0, p)
    d_beta_A_sc = h.beta_A_norm(1, p) - h.beta_A_norm(0, p)

    T = data_dict[p_str] [:, 0]
    Delta2_expt = data_dict[p_str] [:, 6]
    d_Delta2_expt = data_dict[p_str] [:, 7]

    line = ax3.errorbar(T, Delta2_expt, d_Delta2_expt, 
                        c=color_dict[p_str], 
                        ls='', marker='.', label=p_str)

    ax3.plot(t_arr*Tc_mK, (bcs_gap_A_phase_arr*Tc_mK)**2,  
             c=line[0].get_c(), ls=':')
             # label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')

    gap2_hei24 = (-0.5*alpha4/beta4_gap)/tr_AdagA *Tc_mK**2/(sc_corr_factor*hei24_fun(t_arr, p))

    ax3.plot(t_arr*Tc_mK, gap2_hei24, 
             c='k', ls='--')
    
    gap2_this_d_beta_fix = (-0.5*alpha4/(beta4_gap + d_beta_A_sc))/tr_AdagA *Tc_mK**2
    ax3.plot(t_arr*Tc_mK, gap2_this_d_beta_fix, 
             c=line[0].get_c(), ls='--', alpha=0.5)

    d_beta_A_sc = beta_corr_fitfun(t_arr*Tc_mK, *popt_bc_dict[p_str])*h.beta_const*2
    gap2_this = (-0.5*alpha4/(beta_A_bcs_approx(t_arr) + d_beta_A_sc))/tr_AdagA *Tc_mK**2
    ax3.plot(t_arr*Tc_mK, gap2_this, 
             c=line[0].get_c(), ls='--')

ax3.set_title(r'Data (markers), BCS (dot), ' 
              # + ' quartic fit to BCS free energy (dot-dash)'
              # + '\n'
              + r'Heikkinen et al 2024 (black dash), '
              + '\n'
              + r'$\Delta^2 = (1 - t)/2(\beta_A^{\rm BCS}(t) + \Delta\beta_A^{\rm sc})$' 
              + r' ($\Delta\beta_A^{\rm sc}$ fixed: faint dash; fitted: dash)',
              fontsize='smaller')

ax3.set_xlabel(r'$T/T_c$')
ax3.set_ylabel(r'$\Delta^2/k_{\rm B}^2$ [mK$^2$]')

ax3.set_xlim(0,3)
ax3.set_ylim(0,35)
# ax3.set_ylim(0.8,1.6)

ax3.legend(title = r'$p$/bar')
ax3.grid(True)

fig3.tight_layout()

if savefig:
    fig3.savefig(fig_dir + 'gap2_A_data_fits_T.png')
    fig3.savefig(fig_dir + 'gap2_A_data_fits_T.pdf')


