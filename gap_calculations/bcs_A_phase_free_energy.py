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

# import matplotlib_inline;
# matplotlib_inline.backend_inline.set_matplotlib_formats('png', 'jpeg');



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
h.set_default('DEFAULT_SC_CORRS', 'RWS19-interp')

data_dir = '/Users/hindmars/physics/He3/experiments/RHUL_gap/'

file_pattern = 'Run{:}_expr_gap_specular_{:}bar.dat'

run_list = [33]*3 + [34]*2
p_list = [0.17, 2.46, 5.50, 12.00, 21.00]
p_str_list = [f'{p:.2f}' for p in p_list]

color_dict = dict(zip(p_str_list, bcs.color_list))

data_dict = {}

for p_str, run in zip(p_str_list, run_list):
    data_dict[p_str] = np.loadtxt(data_dir + file_pattern.format(run, p_str), skiprows=20)

# #%% Check Tc suppression with Tc from experiment

# Tc_mK_RHUL_dict = {}
# beta_A_RHUL_dict = {}

# fig7, ax7 = plt.subplots()


# def linear(x, *a):
#     return a[0] + a[1]*x

# for p_str in p_str_list: 

#     p = float(p_str)
#     Tc_mK =  h.Tc_mK(p)
    
#     T = data_dict[p_str] [:, 0]
#     Delta2_expt = data_dict[p_str] [:, 6]
#     d_Delta2_expt = data_dict[p_str] [:, 7]

#     just_below_Tc = (T < Tc_mK - 0.02) & (Tc_mK -0.05 < T)
#     near_Tc = (np.abs(T - Tc_mK) < 0.05)
    
#     popt, copt = scopt.curve_fit(linear, T[just_below_Tc], Delta2_expt[just_below_Tc], 
#                                  p0=(-Tc_mK/bcs.bcs_const, 1/bcs.bcs_const))

#     Tc_mK_est = -popt[0]/popt[1]
#     beta_est = -1/popt[1]*Tc_mK_est/4

#     ax7.errorbar(T[near_Tc]- Tc_mK, Delta2_expt[near_Tc], d_Delta2_expt[near_Tc], 
#                  ls='', marker='.', label=p_str + r': $T_c^{{\rm est}} = {:.2f}$ mK, $T_c^{{\rm G}} = {:.2f}$ mK'.format(Tc_mK_est, Tc_mK) )
#     ax7.plot((T[just_below_Tc] - Tc_mK), linear(T[just_below_Tc], *popt), 'k--')

#     print(p_str, Tc_mK, Tc_mK_est, beta_est, h.beta_A_norm(1,p))
    
#     Tc_mK_RHUL_dict[p_str] = Tc_mK_est
#     beta_A_RHUL_dict[p_str] = beta_est

# ax7.set_xlabel(r'$(T - T_c^{\rm G})/$mK')
# ax7.set_ylabel(r'$\Delta^2/(k_{\rm B} T_c)^2$')
# ax7.grid(True)
# ax7.legend(fontsize='smaller')
# ax7.set_xlim(-0.06,0.02)

#%% Plot and fit various approximations to A phase free energy

fig3, ax3 = plt.subplots()

f_lim = 0.2
delta_A_0 = 2.0293

p=5.5

def f_A_fitter(x, t, *b):
    val = (t - 1)*x**2 * 2
    for n, par in enumerate(b):
        val += par*(2*x**2)**(n+2)
    return val

def bcs_gap_A_phase_approx(t, *a): 
    return a[0] * (1 - t**a[1])**0.5

def beta_A_bcs_gap_approx(t, *a): 
    return 0.25*(1 - t)/bcs_gap_A_phase_approx(t, *a)**2

for n in range(15,len(t_arr),8):
    
    f_A_arr = f_A_t_arr[n,:]
    
    t = t_arr[n]
    
    def f_A_fitter_this_t(x, *b):
        return f_A_fitter(x, t, *b)

    
    f_min = np.min(f_A_arr)
    
    in_range = np.logical_and(f_min < f_A_arr, f_A_arr<f_lim)
    
    popt, copt = scopt.curve_fit(f_A_fitter_this_t, (D_arr)[in_range], f_A_arr[in_range], p0=(.1,))

    df_A = f_A_t_arr[n,:] - f_A_t_arr[n,0]    
    line = ax3.plot(D_arr, df_A)

    ax3.plot(D_arr, f_A_fitter_this_t(D_arr, *popt), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.2f}$: '.format(t_arr[n]) + r'$b = {:.3f}$'.format(*popt))

    idx_min = np.argmin(df_A)
    ax3.scatter(D_arr[idx_min], df_A[idx_min], marker='o',  c=line[0].get_c())

    print(bcs_gap_A_phase_arr[n], D_arr[idx_min], *popt)

    v, A_GL_arr, f_GL_arr = h3fe.line_section(h.z3, 'A', t, p, scale=np.max(D_arr))

    D_GL_arr = h.norm(A_GL_arr)
    ax3.plot(D_GL_arr, f_GL_arr, c=line[0].get_c(), ls='-.')

    ax3.scatter(h.delta_A_norm(t, p), h.f_A_norm(t, p), marker='^', c=line[0].get_c())

    ax3.plot(D_GL_arr, bcs.quartic_symmetric(D_GL_arr, (t-1), h.beta_A_norm(0, 0)), 
             c=line[0].get_c(), ls=':' )

    ax3.plot(D_GL_arr, bcs.quartic_symmetric(D_GL_arr*np.sqrt(2), (t-1), beta_A_bcs_gap_approx(t, delta_A_0, 3)), 
             # c='k',
             c=line[0].get_c(), 
             ls=bcs.ls_dict['long dash'] )

ax3.set_xlabel(r'$\Delta/k_B T_c$')
ax3.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax3.set_xlim(0, 5)
ax3.set_ylim(-5, 3)

ax3.legend(fontsize='small')

ax3.set_title(r'BCS A-phase free energy (solid), fit to $2(t-1)x^2/2 + 4bx^4$ (dash), $b=(1-t)/2(\Delta_A^{\rm BCS})^2$ (long dash)),'
              + '\n'
              + 'GL+SC (dot-dash), GL (dot). ' 
              # + '\n'
              + r'$p = {:.1f}$, '.format(p) 
              + r'Fitted $\beta_{{245}} = {:.5f}$  '.format(popt[0])
              + r'$\beta_{{245}}^{{\rm BCS}} = {:.5f}$'.format(h.beta_A_norm(0, 0)), 
              fontsize='smaller')

ax3.grid(True)

#%% Fit A phase free energy to quartic (t-1)gap^2 + beta*gap^4)

# t_arr = np.linspace(0.05,1.15,23)

popt_f_A_list = []


for n, t in enumerate(t_arr):
    
    f_A_arr = f_A_t_arr[n,:]
    
    t = t_arr[n]
    
    f_min = np.min(f_A_arr)
    
    in_range = np.logical_and(f_min < f_A_arr, f_A_arr<f_lim)

    # popt, copt = scopt.curve_fit(quartic_symmetric, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(-t, .1))

    def f_A_fitter_this_t(x, *b):
        return f_A_fitter(x, t, *b)

    popt, copt = scopt.curve_fit(f_A_fitter_this_t, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(.1))
    
    popt_f_A_list.append(popt)

popt_f_A_arr = np.array(popt_f_A_list)

#%% Plot effective quartic parameter 

fig2, ax2 = plt.subplots()

def beta_bcs_approx(t, *a): 
    return 6/(1 + a[0]*(t) + (2 - a[0])*(t)**2)

alpha4 = t_arr - 1
beta4_f_A = popt_f_A_arr[:,0]

beta4_gap = (1-t_arr)/bcs_gap_A_phase_arr**2/4

in_range = (t_arr < 1) & (t_arr > 0.4) 

popt_f_beta_A, copt_f_beta_A = scopt.curve_fit(beta_bcs_approx, 
                                           t_arr[in_range], beta4_f_A[in_range]/h.beta_const, 
                                           p0=(1,))

popt_gap_beta_A, copt_gap_beta_A = scopt.curve_fit(beta_bcs_approx, 
                                           t_arr[in_range], beta4_gap[in_range]/h.beta_const, 
                                           p0=(1,))

line_beta_f = ax2.plot(t_arr, beta4_f_A/h.beta_const, 
                       c='b',
                       label=r'Fit to $2(t - 1)\Delta^2 + 4\beta\Delta^4$')
                       # label=r'$\tilde\beta(t)/\beta_0$')

ax2.plot(t_arr, beta_bcs_approx(t_arr, *popt_f_beta_A), 
         ls='--', c=line_beta_f[0].get_c(), 
         label=r'$6/(1 + a_0t + (2-a_0)t^2)$' + '\n' + r'$a_0 = {:.2f}$'.format(*popt_f_beta_A))

line_beta_gap = ax2.plot(t_arr[t_arr<1], beta4_gap[t_arr<1]/h.beta_const, 
         # ls=bcs.ls_dict['long dash'], 
         c='r', 
         # label=r'$(1-t)/4\Delta^2_{\rm BCS}$')
         label=r'$\beta_A^{\rm BCS} =  = (1 - t)/2 [\Delta^{BCS}_A(t)]^2$')

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

    # Tc_mK = Tc_mK_RHUL_dict[p_str]
    Tc_mK = h.Tc_mK(p)
    
    t_expt = T/Tc_mK
    beta4_expt = (1-t_expt)/(Delta2_expt/Tc_mK**2)/4
    d_beta4_expt = (1-t_expt)/((Delta2_expt - d_Delta2_expt)/Tc_mK**2)/4 - beta4_expt
    
    
    # line_beta_expt = ax2.errorbar(t_expt[t_expt<1.0], beta4_expt[t_expt<1.0]/h.beta_const,
    #                              d_beta4_expt[t_expt<1.0]/h.beta_const, 
    #                              marker='.',
    #                              markersize=1,
    #                              ls='',
    #                              c=color_dict[p_str],
    #          label=r'$p$/bar = ' + p_str)
    # ax2.scatter(1, h.beta_A_norm(1,p)/h.beta_const, 
    #             c=line_beta_expt[0].get_c())


ax2.set_title(r'Effective quartic parameter $\beta_A^{\rm eff}$')
              # + '\n'
              # + 'Expt (errorbar), RWS19($T_c$) (dot)')

ax2.set_xlabel(r'$T/T_c$')
ax2.set_ylabel(r'$\beta_A^{\rm eff}(t)/\beta^{\rm wc}$')

ax2.set_xlim(0,1.1)
ax2.set_ylim(1.0, 4.0)

ax2.legend(fontsize='small')
ax2.grid(True)


#%% Fit free energy to sextic


fig3, ax3 = plt.subplots()

for n in range(7,len(t_arr),8):
    f_A_arr = f_A_t_arr[n,:]
    
    t = t_arr[n]
    
    # def f_A_fitter_this_t(x, *b):
    #     return f_A_fitter(x, t, *b)

    this_fit_fun = bcs.poly_in_squares

    f_min = np.min(f_A_arr)
    
    in_range = np.logical_and(f_min < f_A_arr, f_A_arr<f_lim)
    
    popt6, copt6 = scopt.curve_fit(this_fit_fun, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(t-1, .1, -.01))
    
    line = ax3.plot(D_arr, f_A_arr - f_A_arr[0])

    ax3.plot(D_arr, this_fit_fun(D_arr, *popt6), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.1f}$: '.format(t_arr[n]) + r'$a = ({:.2f}, {:.2f})$'.format(*popt6))


ax3.set_xlabel(r'$\Delta/k_B T_c$')
ax3.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax3.set_xlim(0, 3)
ax3.set_ylim(-5, 5)

ax3.legend(fontsize='small')

ax3.set_title('BCS A-phase free energy (solid), fit to $a_0x^2 + a_1x^4 + a_2x^6$ (dash)')

ax3.grid(True)


#%% Fit free energy to sextic

popt6_list = []
copt6_list = []

for n, t in enumerate(t_arr):
    
    f_A_arr = f_A_t_arr[n,:]
    
    t = t_arr[n]
    
    # def f_A_fitter_this_t(x, *b):
    #     return f_A_fitter(x, t, *b)

    this_fit_fun = bcs.poly_in_squares

    f_min = np.min(f_A_arr)
    
    in_range = np.logical_and(f_min < f_A_arr, f_A_arr<f_lim)
    
    popt6, copt6 = scopt.curve_fit(this_fit_fun, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(t-1, .1, -.01))
    
    popt6_list.append(popt6)
    copt6_list.append(copt6)

popt6_arr = np.array(popt6_list)
copt6_arr = np.array(copt6_list)
#%%

fig4, ax4 = plt.subplots()

alpha6 = popt6_arr[:,0]/2
beta6 = popt6_arr[:,1]/4
gamma6 = popt6_arr[:,2]/8

beta4_f_A = popt_f_A_arr[:,0]

ind_Tc = np.argmin(abs(t_arr - 1))
gamma6_Tc = gamma6[ind_Tc]

ax4.plot(t_arr[:-1], -(alpha6)[:-1], label=r'$-\tilde\alpha(t)$')

ax4.plot(t_arr, beta4_f_A/h.beta_const, label=r'$\tilde\beta(t)/\beta^{\rm wc}$ (quartic fit)')

line_beta = ax4.plot(t_arr, beta6/h.beta_const, label=r'$\tilde\beta(t)/\beta^{\rm wc}$')
# ax4.plot(t_arr, h.beta_A_norm(np.zeros_like(t_arr),0)/h.beta_const/t_arr**2, ls='--', c=line_beta[0].get_c(), label=r'$\beta_{245}/\beta_0 t^2$')
# ax4.plot(t_arr, 2/t_arr**2, ls='--')

line_gamma6 = ax4.plot(t_arr, gamma6/h.beta_const, label=r'$-\tilde\gamma(t)/(\beta^{{\rm wc}})$, $\gamma(T_c) = {:.3g}$'.format(gamma6_Tc))

ax4.set_title(r'BCS free energy fit parameters: ' + r'$f = (t - 1)\Delta^2 + \beta\Delta^4+ \gamma\Delta^6$')


ax4.set_xlabel(r'$T/T_c$')

ax4.set_xlim(0,1.25)
ax4.set_ylim(-1, 7)

ax4.legend()
ax4.grid(True)

# #%%

# fig5, ax5 = plt.subplots()

# approx_gap2_A_plus, approx_gap2_A_minus  = bcs.solve_quadratic(3*gamma6, 2*beta6, alpha6)

# p = 22
# # delta_b_B = np.sum(h.delta_b_asarray(p, [1,2])) + np.sum(h.delta_b_asarray(p, [3,4,5]))/3 
# # delta_beta_B = h.beta_const * delta_b_B
# delta_b_A = np.sum(h.delta_b_asarray(p, [2, 4, 5])) 
# delta_beta_A = h.beta_const * delta_b_A

# sc_corr_factor = h.beta_A_norm(t_arr, p)/h.beta_A_norm(0, p)
# approx_gap2_A_plus_sc, _  = bcs.solve_quadratic(3*gamma6*sc_corr_factor, 2*beta6*sc_corr_factor, alpha6)

# # approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6_Tc*np.ones_like(t_arr)/t_arr**3, 2*h.beta_A_norm(0*t_arr, 0)/t_arr**2, 
# #                                                        t_arr - 1 - 0.5*(t_arr - 1)**2)


# ax5.plot(t_arr, approx_gap2_A_plus/2, 
#          label=r'$\Delta_+(\alpha,\beta,\gamma)/k_BT_c$')
# ax5.plot(t_arr, approx_gap2_A_plus_sc/2, 
#          label=r'$\Delta_+(\alpha,\beta^{\rm sc},\gamma)/k_BT_c$')
# # ax5.plot(t_arr, np.sqrt(approx_gap2_minus)/np.sqrt(3), label=r'$\Delta_-$')


# ax5.plot(t_arr, bcs_gap_A_phase_arr**2, 'b--', label=r'$\Delta_{\rm A, BCS}(t)/k_BT_c$')


# ax5.set_xlabel(r'$T/T_c$')
# ax5.set_ylabel(r'$(\Delta/k_BT_c)^2$')

# ax5.set_xlim(0,1.25)
# ax5.set_ylim(0,5)

# ax5.legend()
# ax5.grid(True)


#%% Fit to even poly up to $D^8$

fig7, ax7 = plt.subplots()

for n in range(7,len(t_arr),8):
    
    f_A_arr = f_A_t_arr[n,:]
    
    t = t_arr[n]
    
    f_min = np.min(f_A_arr)
    
    in_range = np.logical_and(f_min < f_A_arr, f_A_arr<f_lim)
    
    popt, copt = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(-t, .1, .1, .1))
    
    line = ax7.plot(D_arr, f_A_arr - f_A_arr[0])

    ax7.plot(D_arr, bcs.poly_in_squares(D_arr, *popt), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.1f}$: '.format(t_arr[n]) + r'$a = ({:.2f}, {:.2f}, {:.3f})$'.format(*popt))


ax7.set_xlabel(r'$\Delta/k_B T_c$')
ax7.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax7.set_xlim(0, 3)
ax7.set_ylim(-5, 5)

ax7.legend(fontsize='small')

ax7.set_title('BCS A-phase free energy (solid), fit to $a_0x^2 + a_1x^4 + a_2x^6 + a_3x^8$ (dash)')

ax7.grid(True)

#%% extract fitting parameters for $D^8$

popt8_list = []
copt8_list = []

for n, t in enumerate(t_arr):
    
    f_A_arr = f_A_t_arr[n,:]
    
    # t = t_arr[n]
    
    popt8, copt8 = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_A_arr<f_lim], f_A_arr[f_A_arr<f_lim], p0=(-t, .1, -.1, .1))
    
    popt8_list.append(popt8)
    copt8_list.append(copt8)

popt8_arr = np.array(popt8_list)
copt8_arr = np.array(copt8_list)

#%% plot fitting parameters for $D^8$

fig8, ax8 = plt.subplots()

alpha8 = popt8_arr[:,0]/2
beta8 = popt8_arr[:,1]/4
gamma8 = popt8_arr[:,2]/8
delta8 = popt8_arr[:,3]/16

ind_Tc = np.argmin(abs(t_arr - 1))
gamma8_Tc = gamma8[ind_Tc]

ax8.plot(t_arr[:-1], -(alpha6)[:-1], label=r'$-\tilde\alpha(t)$')

line_beta8 = ax8.plot(t_arr, beta8/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')
ax8.plot(t_arr, h.beta_A_norm(np.zeros_like(t_arr),0)/h.beta_const/t_arr**2, 
         ls='--', c=line_beta8[0].get_c(), label=r'$\beta_{245}/\beta_0 t^2$')

line_gamma8 = ax8.plot(t_arr, gamma8/gamma8_Tc, label=r'$-\tilde\gamma(t)/\gamma(T_c)$')
# ax8.plot(t_arr, 1/t_arr**3, 
#          ls='--', c=line_gamma8[0].get_c(), label=r'$\gamma(T_c) / t^3$')


ax8.set_title(r'BCS free energy fit parameters: ' + r'$f = \alpha\Delta^2 + \beta\Delta^4+ \gamma\Delta^6+ \delta\Delta^8$')


ax8.set_xlabel(r'$T/T_c$')

ax8.set_xlim(0,1.25)
ax8.set_ylim(0,10)

ax8.legend()
ax8.grid(True)

#%% plot gap against temperature for D^8 model

from sympy.solvers import solve
from sympy import symbols

x = symbols('x')

# a, b, c, d = symbols('alpha beta gamma delta')

def solve_cubic_arr(*a_arrs):
    
    sols = []
    for a in zip(*a_arrs):
        sols.append(solve(a[0] + a[1]*x + a[2]*x**2 + a[3]*x**3, x))
        
    return np.array(sols)


fig9, ax9 = plt.subplots()

# approx_gap2_A_plus, approx_gap2_A_minus  = solve_quadratic(3*gamma8, 2*beta8, alpha8)
# approx_gap2_A_plus  = solve_quadratic(3*gamma8, 2*beta8, alpha8)[0]/2
approx_gap2_A_plus  = solve_cubic_arr(alpha8, 2*beta8, 3*gamma8, 4*delta8)[:,0]/2

p = 22
# delta_b_B = np.sum(h.delta_b_asarray(p, [1,2])) + np.sum(h.delta_b_asarray(p, [3,4,5]))/3 
# delta_beta_B = h.beta_const * delta_b_B
delta_b_A = np.sum(h.delta_b_asarray(p, [2, 4, 5])) 
delta_beta_A = h.beta_const * delta_b_A

# approx_gap2_A_plus_sc  = solve_quadratic(3*gamma8, 2*(beta8 + t_arr*delta_beta_A), alpha8)[0]/2
approx_gap2_A_plus_sc  = solve_cubic_arr(alpha8, 2*(beta8+ t_arr*delta_beta_A), 3*gamma8, 4*delta8)[:,0]/2

# approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6_Tc*np.ones_like(t_arr)/t_arr**3, 2*h.beta_A_norm(0*t_arr, 0)/t_arr**2, 
#                                                        t_arr - 1 - 0.5*(t_arr - 1)**2)


ax9.plot(t_arr, approx_gap2_A_plus, 
         label=r'$\Delta_+(\alpha,\beta,\gamma)/k_BT_c$')
ax9.plot(t_arr, approx_gap2_A_plus_sc, 
         label=r'$\Delta_+(\alpha,\beta^{\rm sc},\gamma)/k_BT_c$')
# ax5.plot(t_arr, np.sqrt(approx_gap2_minus)/np.sqrt(3), label=r'$\Delta_-$')

ax9.plot(t_arr, bcs_gap_A_phase_arr**2, 'b--', label=r'$\Delta_{\rm A, BCS}(t)/k_BT_c$')


ax9.set_xlabel(r'$T/T_c$')
ax9.set_ylabel(r'$(\Delta/k_BT_c)^2$')

ax9.set_xlim(0,1.25)
ax9.set_ylim(0,5)

ax9.legend()
ax9.grid(True)

#%% 

# fig10, ax10 = plt.subplots()

# foo = h3fe.line_section(h.z3, 'A', 0.5, 5.5, scale=np.max(D_arr))

# def plot_GL_potential()

