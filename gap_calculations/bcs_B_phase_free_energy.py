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


bcs_gap_B_dict = np.load('bcs_gap_B_phase.npz')

bcs_f_B_dict = np.load('free_energy_B_phase.npz')

t_arr = bcs_f_B_dict['t']
D_arr = bcs_f_B_dict['delta_B']
f_B_t_arr = bcs_f_B_dict['free_energy_B_phase']

bcs_gap_B_phase_arr = bcs_gap_B_dict['bcs_gap_B_phase']

#%% Plot and fit various approximations to b phase free energy

def f_B_fitter(x, t, *b):
    val = (t - 1)*x**2 * 2
    for n, par in enumerate(b):
        val += par*(2*x**2)**(n+2)
    return val


fig3, ax3 = plt.subplots()

f_lim = 0.05

p=5.5

for n in range(15,len(t_arr),8):
    
    f_B_arr = f_B_t_arr[n,:]
    
    t = t_arr[n]
    
    # def f_B_fitter_this_t(x, *b):
    #     return f_B_fitter(x, b[0], *b[1:])

    fit_fun = bcs.poly_in_squares

    f_min = np.min(f_B_arr)
    
    in_range = np.logical_and(f_min < f_B_arr, f_B_arr<f_lim)
    
    popt, copt = scopt.curve_fit(fit_fun, (D_arr)[in_range], f_B_arr[in_range], p0=(t, .1,))

    df_B = f_B_t_arr[n,:] - f_B_t_arr[n,0]    
    line = ax3.plot(D_arr, df_B)

    ax3.plot(D_arr, fit_fun(D_arr, *popt), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.2f}$: '.format(t_arr[n]) + r'$a = ({:.3f}, {:.3f})$'.format(*popt) )
    idx_min = np.argmin(df_B)
    ax3.scatter(D_arr[idx_min], df_B[idx_min], marker='o',  c=line[0].get_c())

    print(bcs_gap_B_phase_arr[n], D_arr[idx_min]) #marker='s', c=line[0].get_c())    

    foo = h3fe.line_section(h.z3, 'B', t, p, scale=np.max(D_arr))

    ax3.plot(h.norm(foo[1]), foo[2], c=line[0].get_c(), ls=':')

    ax3.scatter(h.delta_B_norm(t, p), h.f_B_norm(t, p), marker='^', c=line[0].get_c())

    # print(h.f_B_norm(t, 0), np.min(foo[2]))
    # print(popt[0]/2, popt[1]/4)

ax3.set_xlabel(r'$\Delta/k_B T_c$')
ax3.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax3.set_xlim(0, 4)
ax3.set_ylim(-4, 4)

ax3.legend(fontsize='small')

ax3.set_title('BCS B-phase free energy (solid), fit to $a_0x^2 + a_1x^4$ (dash), GL (dot)' 
              + '\n'
              + r'Fitted $\beta_{{B}} = {:.5f}$  '.format(popt[1]/9)
              + r'$\beta_{{B}}^{{\rm BCS}} = {:.5f}$'.format(h.beta_B_norm(0, 0)), 
              fontsize='smaller')

ax3.grid(True)

#%% Tit quartic approximations to b phase free energy


# t_arr = np.linspace(0.05,1.15,23)

popt_list = []

for n, t in enumerate(t_arr):
    
    f_B_arr = f_B_t_arr[n,:]
    
    t = t_arr[n]
    
    # def f_B_fitter_this_t(x, *b):
    #     return f_B_fitter(x, b[0], *b[1:])
    
    
    fit_fun = bcs.poly_in_squares

    f_min = np.min(f_B_arr)
    
    in_range = np.logical_and(f_min < f_B_arr, f_B_arr<f_lim)

    popt, copt = scopt.curve_fit(fit_fun, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(.1,.1,))
    # popt, copt = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(-t, .1))
    
    popt_list.append(popt)

popt_arr = np.array(popt_list)

#%% Plot effective beta parameter


def bcs_gap_B_phase_approx(t, *a): 
    return a[0] * (1 - t**a[1])**0.5

def beta_B_bcs_gap_approx(t, *a): 
    return 0.25*(1 - t)/bcs_gap_B_phase_approx(t, *a)**2

def beta_B_bcs_approx(t, *a):
    return a[0] + a[1]*(t-1)


fig2, ax2 = plt.subplots()

alpha_f_B = popt_arr[:,0]/3
# alpha = - 1 + t_arr
beta_f_B = popt_arr[:,1]/9

line_alpha = ax2.plot(t_arr, -alpha_f_B)

line_beta = ax2.plot(t_arr, beta_f_B/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')

ax2.plot(t_arr, h.beta_B_norm(np.zeros_like(t_arr),0)/h.beta_const, 
         ls='--', c='k',
         label=r'$\tilde\beta/\beta_0 = 5/3$')

ax2.plot(t_arr, 1 - t_arr, ls='--', c = line_alpha[0].get_c(),
         label=r'$1- t$')

popt_f_beta_B, copt_f_beta_B = scopt.curve_fit(beta_B_bcs_approx, 
                                           t_arr[t_arr>0.5], beta_f_B[t_arr>0.5]/h.beta_const, 
                                           p0=(1,1))

ax2.plot(t_arr, beta_B_bcs_approx(t_arr, *popt_f_beta_B), 
         ls='--', c=line_beta[0].get_c(),
         label='Fit: ' +  r'$a = ({:.3f}, {:.3f})$'.format(*popt_f_beta_B))

# ax2.plot(t_arr, np.log(t_arr)/(t_arr-1), ls='--', c='b', label=r'$\ln(t)/(t -1)$')

# ax2.plot(t_arr, bcs_gap_arr, 'b--', label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')

ax2.set_title(r'BCS free energy fit parameters: ' + r'$f = \alpha\Delta^2 + \beta\Delta^4$')

ax2.set_xlabel(r'$T/T_c$')

ax2.set_xlim(0,1.25)
ax2.set_ylim(0,5)

ax2.legend()
ax2.grid(True)

#%%

popt6_list = []
copt6_list = []

for n, t in enumerate(t_arr):
    
    f_B_arr = f_B_t_arr[n,:]
    
    t = t_arr[n]
    
    f_min = np.min(f_B_arr)
    
    in_range = np.logical_and(f_min < f_B_arr, f_B_arr<f_lim)
    
    popt6, copt6 = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(-t, .1, -.01))
    
    popt6_list.append(popt6)
    copt6_list.append(copt6)

popt6_arr = np.array(popt6_list)
copt6_arr = np.array(copt6_list)
#%% Sextic 

fig4, ax4 = plt.subplots()

alpha6 = popt6_arr[:,0]/3
beta6 = popt6_arr[:,1]/9
gamma6 = popt6_arr[:,2]/27

ind_Tc = np.argmin(abs(t_arr - 1))
gamma6_Tc = gamma6[ind_Tc]

ax4.plot(t_arr[:-1], -(alpha6)[:-1], label=r'$-\tilde\alpha(t)$')
# ax2.plot(t_arr, (t_arr - 1), 'k--')

line_beta = ax4.plot(t_arr, beta6/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')
ax4.plot(t_arr, h.beta_B_norm(np.zeros_like(t_arr),0)/h.beta_const/t_arr**2, ls='--', c=line_beta[0].get_c(), label=r'$\beta_{B}/\beta_0 t^2$')
# ax4.plot(t_arr, 0.8/t_arr**2, ls='--')
# ax4.plot(t_arr, np.log(t_arr)/(t_arr-1), ls='--', c='b')

line_gamma6 = ax4.plot(t_arr, gamma6/gamma6_Tc, label=r'$-\tilde\gamma(t)/\gamma(T_c)$')

ax4.set_title(r'BCS free energy fit parameters: ' + r'$f = \alpha\Delta^2 + \beta\Delta^4+ \gamma\Delta^6$')


# ax4.plot(t_arr, np.sqrt(-0.5*alpha6/beta6 ), 'k', label=r'$\sqrt{-\tilde{\alpha}/2\tilde{\beta}}$')
# ax4.plot(t_arr, bcs_gap_arr, 'b--', label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')

ax4.set_xlabel(r'$T/T_c$')

ax4.set_xlim(0,1.25)
ax4.set_ylim(0,10)

ax4.legend()
ax4.grid(True)

#%%

fig5, ax5 = plt.subplots()

approx_gap2_B_plus, approx_gap2_B_minus  = bcs.solve_quadratic(3*gamma6, 2*beta6, alpha6)

p = 22
delta_b_B = np.sum(h.delta_b_asarray(p, [1,2])) + np.sum(h.delta_b_asarray(p, [3,4,5]))/3 
delta_beta_B = h.beta_const * delta_b_B
# delta_b_B = np.sum(h.delta_b_asarray(p, [3, 4, 5])) 
# delta_beta_B = h.beta_const * delta_b_B

sc_corr_factor = h.beta_B_norm(t_arr, p)/h.beta_B_norm(0, p)
approx_gap2_B_plus_sc, _  = bcs.solve_quadratic(3*gamma6*sc_corr_factor, 2*beta6*sc_corr_factor, alpha6)

# approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6_Tc*np.ones_like(t_arr)/t_arr**3, 2*h.beta_B_norm(0*t_arr, 0)/t_arr**2, 
#                                                        t_arr - 1 - 0.5*(t_arr - 1)**2)


ax5.plot(t_arr, approx_gap2_B_plus/2, 
         label=r'$\Delta_+(\alpha,\beta,\gamma)/k_BT_c$')
ax5.plot(t_arr, approx_gap2_B_plus_sc/2, 
         label=r'$\Delta_+(\alpha,\beta^{\rm sc},\gamma)/k_BT_c$')
# ax5.plot(t_arr, np.sqrt(approx_gap2_minus)/np.sqrt(3), label=r'$\Delta_-$')


ax5.plot(t_arr, bcs_gap_B_phase_arr**2, 'b--', label=r'$\Delta_{\rm B, BCS}(t)/k_BT_c$')


ax5.set_xlabel(r'$T/T_c$')
ax5.set_ylabel(r'$(\Delta/k_BT_c)^2$')

ax5.set_xlim(0,1.25)
ax5.set_ylim(0,5)

ax5.set_title('B-phase gap: 6th order fit comparisons (with SC model)')


ax5.legend()
ax5.grid(True)


    
#%% Fit to even poly up to $D^8$

fig7, ax7 = plt.subplots()

for n in range(7,len(t_arr),8):
    
    f_B_arr = f_B_t_arr[n,:]
    
    t = t_arr[n]
    
    f_min = np.min(f_B_arr)
    
    in_range = np.logical_and(f_min < f_B_arr, f_B_arr<f_lim)
    
    popt, copt = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(-t, .1, .1, .1))
    
    line = ax7.plot(D_arr, f_B_arr - f_B_arr[0])

    ax7.plot(D_arr, bcs.poly_in_squares(D_arr, *popt), 
            c=line[0].get_c(), ls='--', 
            label=r'$t = {:.1f}$: '.format(t_arr[n]) + r'$a = ({:.2f}, {:.2f}, {:.3f})$'.format(*popt))


ax7.set_xlabel(r'$\Delta/k_B T_c$')
ax7.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax7.set_xlim(0, 3)
ax7.set_ylim(-5, 5)

ax7.legend(fontsize='small')

ax7.set_title('BCS B-phase free energy (solid), fit to $a_0x^2 + a_1x^4 + a_2x^6 + a_3x^8$ (dash)')

ax7.grid(True)

#%% extract fitting parameters for $D^8$

popt8_list = []
copt8_list = []

for n, t in enumerate(t_arr):
    
    f_B_arr = f_B_t_arr[n,:]
    
    # t = t_arr[n]
    
    popt8, copt8 = scopt.curve_fit(bcs.poly_in_squares, (D_arr)[f_B_arr<f_lim], f_B_arr[f_B_arr<f_lim], p0=(-t, .1, -.1, .1))
    
    popt8_list.append(popt8)
    copt8_list.append(copt8)

popt8_arr = np.array(popt8_list)
copt8_arr = np.array(copt8_list)

#%% plot fitting parameters for $D^8$

fig8, ax8 = plt.subplots()

alpha8 = popt8_arr[:,0]/3
beta8 = popt8_arr[:,1]/3**2
gamma8 = popt8_arr[:,2]/3**3
delta8 = popt8_arr[:,3]/3**4

ind_Tc = np.argmin(abs(t_arr - 1))
gamma8_Tc = gamma8[ind_Tc]

ax8.plot(t_arr[:-1], -(alpha6)[:-1], label=r'$-\tilde\alpha(t)$')
# ax2.plot(t_arr, (t_arr - 1), 'k--')

line_beta8 = ax8.plot(t_arr, beta8/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')
ax8.plot(t_arr, h.beta_B_norm(np.zeros_like(t_arr),0)/h.beta_const/t_arr**2, 
         ls='--', c=line_beta8[0].get_c(), label=r'$\beta_{B}/\beta_0 t^2$')

line_gamma8 = ax8.plot(t_arr, gamma8/gamma8_Tc, label=r'$-\tilde\gamma(t)/\gamma(T_c)$')
ax8.plot(t_arr, 1/t_arr**3, 
         ls='--', c=line_gamma8[0].get_c(), label=r'$\gamma(T_c) / t^3$')


ax8.set_title(r'BCS free energy fit parameters: ' + r'$f = \alpha\Delta^2 + \beta\Delta^4+ \gamma\Delta^6+ \delta\Delta^8$')


ax8.set_xlabel(r'$T/T_c$')

ax8.set_xlim(0,1.25)
ax8.set_ylim(0,10)

ax8.legend()
ax8.grid(True)

#%% plot gap against temperature for D^8 model


fig9, ax9 = plt.subplots()

# approx_gap2_B_plus, approx_gap2_B_minus  = solve_quadratic(3*gamma8, 2*beta8, alpha8)
# approx_gap2_B_plus  = solve_quadratic(3*gamma8, 2*beta8, alpha8)[0]/2
approx_gap2_B_plus  = bcs.solve_cubic_arr(alpha8, 2*beta8, 3*gamma8, 4*delta8)[:,0]/2

p = 22
# delta_b_B = np.sum(h.delta_b_asarray(p, [1,2])) + np.sum(h.delta_b_asarray(p, [3,4,5]))/3 
# delta_beta_B = h.beta_const * delta_b_B
delta_b_B = np.sum(h.delta_b_asarray(p, [2, 4, 5])) 
delta_beta_B = h.beta_const * delta_b_B

# approx_gap2_B_plus_sc  = solve_quadratic(3*gamma8, 2*(beta8 + t_arr*delta_beta_B), alpha8)[0]/2
approx_gap2_B_plus_sc  = bcs.solve_cubic_arr(alpha8, 2*(beta8+ t_arr*delta_beta_B), 3*gamma8, 4*delta8)[:,0]/2

# approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6_Tc*np.ones_like(t_arr)/t_arr**3, 2*h.beta_B_norm(0*t_arr, 0)/t_arr**2, 
#                                                        t_arr - 1 - 0.5*(t_arr - 1)**2)


ax9.plot(t_arr, approx_gap2_B_plus, 
         label=r'$\Delta_+(\alpha,\beta,\gamma)/k_BT_c$')
ax9.plot(t_arr, approx_gap2_B_plus_sc, 
         label=r'$\Delta_+(\alpha,\beta^{\rm sc},\gamma)/k_BT_c$')
# ax5.plot(t_arr, np.sqrt(approx_gap2_minus)/np.sqrt(3), label=r'$\Delta_-$')

ax9.plot(t_arr, bcs_gap_B_phase_arr**2, 'b--', label=r'$\Delta_{\rm B, BCS}(t)/k_BT_c$')


ax9.set_xlabel(r'$T/T_c$')
ax9.set_ylabel(r'$(\Delta/k_BT_c)^2$')

ax9.set_xlim(0,1.25)
ax9.set_ylim(0,5)

ax9.set_title('B-phase gap: 8th order fit comparisons (with SC model)')

ax9.legend()
ax9.grid(True)

#%% 

# fig10, ax10 = plt.subplots()

# foo = h3fe.line_section(h.z3, 'A', 0.5, 5.5, scale=np.max(D_arr))

# def plot_GL_potential()

