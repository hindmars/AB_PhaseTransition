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

xmax = np.inf # Cut-off on integrations
xmin = 0.0
const = np.sqrt(8*np.pi**2/(7*scspe.zeta(3)))

def integrand_1(x, y, t=1):
    e = np.sqrt(x**2 + y**2)
    return 0.5*np.tanh(0.5*e/t)/e

def full_integrand_1(x, y, t=1):
    return integrand_1(x, 0, t) - integrand_1(x, y, t) 

def integrand_2(x, y, t):
    e = np.sqrt(x**2 + y**2)
    return - 0.25 * x**2 * np.cosh(0.5*e/t)**(-2)

def full_integrand_2(x, y, t=1):
    return integrand_2(x, 0, t) - integrand_2(x, y, t) 

def free_energy(D, t):
    """
    D is normlaised eigenvalue of gap matrix
    
    Single D: eigenvalues of the gap matrix assumed equal
    
    Future
    
    - Pair of Ds, evaluate each and add
    - ? Give vector d and compute eigenvalues. Take angular average?
    
    
    """
    
    integral_1 = 2*scint.quad(full_integrand_1, xmin, xmax, args=(D, t))[0]
    integral_2 = 2*scint.quad(full_integrand_2, xmin, xmax, args=(D, t))[0]
    
    return 3*((integral_1 + np.log(t))*D**2 - (integral_2/t + 0.5*D**2) )

def quartic_symmetric(x, *a):
    return a[0]*x**2 + a[1]*x**4

def sextic_symmetric(x, *a):
    return a[0]*x**2 + a[1]*x**4 + a[2]*x**6


# def full_integrand(x, D, t):
#     return integrand_1(x, D/t) - integrand_1(x, 0) 

def gap_eqn(D, t):
    def fun(x):
        # return - integrand_1(x, 0, 1) + integrand_1(x, D, t) 
        return 2*full_integrand_1(x, D, t) 
    res, _ = scint.quad(fun, xmin, xmax)  
    return res + np.log(t)

def gap_GL(t): 
    return const * np.sqrt(1-t)


def free_energy_A_phase(D, t):
    
    mu = np.linspace(-1,1,100)
    gap_arr = D * np.sqrt(1 - mu**2)
    
    f_arr = np.array([free_energy(gap, t) for gap in gap_arr])/2
    
    return np.trapz(f_arr, mu)




#%%

fig, ax = plt.subplots()

# t_list = [0.0625, 0.125, 0.25, 0.5, 0.75, 1.0, 1.25]
t_list = [0.5, 0.75, 1.0, 1.25]
f_lim = 0.1

D_arr = np.linspace(0.01, 3, 100, endpoint=False)

for t in t_list:
    
    f_arr = np.zeros_like(D_arr)
    
    for n, D in enumerate(D_arr):
        f_arr[n] = free_energy(D, t)
    
    popt, copt = scopt.curve_fit(quartic_symmetric, (D_arr)[f_arr<f_lim], f_arr[f_arr<f_lim], p0=(-t, .1))
    
    line = ax.plot(D_arr, f_arr - f_arr[0], label=str(t))

    ax.plot(D_arr, quartic_symmetric(D_arr, *popt), 
            c=line[0].get_c(), ls='--', label=r'$a = ({:.2f}, {:.2f})$'.format(*popt))


ax.set_xlabel(r'$\Delta/k_B T_c$')
ax.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax.set_xlim(0, 5)
ax.set_ylim(-5, 5)

ax.legend(title=r'$T/T_c$', fontsize='small')

ax.set_title('BCS free energy (solid), fit to $a_0x^2 + a_1x^4$ (dash)')

ax.grid(True)

#%%

t_arr = np.arange(0.2,1.2,0.01)

popt_list = []

D_arr = np.linspace(0.01, 2, 100, endpoint=False)

for t in t_arr:
    
    f_arr = np.zeros_like(D_arr)
    
    for n, D in enumerate(D_arr):
        f_arr[n] = free_energy(D, t)
    
    popt, copt = scopt.curve_fit(quartic_symmetric, (D_arr)[f_arr<f_lim], f_arr[f_arr<f_lim], p0=(-t, .1))
    
    popt_list.append(popt)

popt_arr = np.array(popt_list)

#%%
t_arr = np.arange(0.2,1.2,0.01)

bcs_gap_arr = np.array([float(scopt.fsolve(gap_eqn, 1, args=(t,))) for t in t_arr])

#%%

fig2, ax2 = plt.subplots()



alpha = popt_arr[:,0]/3
beta = popt_arr[:,1]/9

ax2.plot(t_arr[:-1], (alpha/(t_arr - 1))[:-1], label=r'$\tilde\alpha(t)/(t -1)$')
# ax2.plot(t_arr, (t_arr - 1), 'k--')

line_beta = ax2.plot(t_arr, beta/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')
ax2.plot(t_arr, h.beta_B_norm(np.zeros_like(t_arr),0)/h.beta_const, ls='--', c=line_beta[0].get_c())

ax2.plot(t_arr, (5/3)/t_arr, ls='--')
ax2.plot(t_arr, np.log(t_arr)/(t_arr-1), ls='--', c='b')


# ax2.plot(t_arr, np.sqrt(-0.5*alpha/beta ), 'k', label=r'$\sqrt{-\tilde{\alpha}/2\tilde{\beta}}$')
# ax2.plot(t_arr, bcs_gap_arr, 'b--', label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')


ax2.set_xlabel(r'$T/T_c$')

ax2.set_xlim(0,1.25)
ax2.set_ylim(0,5)

ax2.legend()
ax2.grid(True)

#%%

fig3, ax3 = plt.subplots()

# t_list = [0.0625, 0.125, 0.25, 0.5, 0.75, 1.0, 1.25]
t_list = [0.5, 0.75, 1.0, 1.25]
f_lim = 0.1

D_arr = np.linspace(0.01, 3, 100, endpoint=False)

for t in t_list:
    
    f_arr = np.zeros_like(D_arr)
    
    for n, D in enumerate(D_arr):
        f_arr[n] = free_energy(D, t)
    
    popt, copt = scopt.curve_fit(sextic_symmetric, (D_arr)[f_arr<f_lim], f_arr[f_arr<f_lim], p0=(-t, .1, .1))
    
    line = ax3.plot(D_arr, f_arr - f_arr[0], label=str(t))

    ax3.plot(D_arr, sextic_symmetric(D_arr, *popt), 
            c=line[0].get_c(), ls='--', label=r'$a = ({:.2f}, {:.2f}, {:.2f})$'.format(*popt))


ax3.set_xlabel(r'$\Delta/k_B T_c$')
ax3.set_ylabel(r'$\Delta f/[\frac{1}{3}N(0)(k_B T_c)^2]$')

ax3.set_xlim(0, 5)
ax3.set_ylim(-5, 5)

ax3.legend(title=r'$T/T_c$', fontsize='small')

ax3.set_title('BCS free energy (solid), fit to $a_0x^2 + a_1x^4 + a_2x^6$ (dash)')

ax3.grid(True)

#%%

t_arr = np.arange(0.2,1.2,0.01)

popt6_list = []
copt6_list = []

D_arr = np.linspace(0.01, 2, 100, endpoint=False)

for t in t_arr:
    
    f_arr = np.zeros_like(D_arr)
    
    for n, D in enumerate(D_arr):
        f_arr[n] = free_energy(D, t)
    
    popt6, copt6 = scopt.curve_fit(sextic_symmetric, (D_arr)[f_arr<f_lim], f_arr[f_arr<f_lim], p0=(-t, .1, -.01))
    
    popt6_list.append(popt6)
    copt6_list.append(copt6)

popt6_arr = np.array(popt6_list)
copt6_arr = np.array(copt6_list)
#%%

fig4, ax4 = plt.subplots()

alpha6 = popt6_arr[:,0]/3
beta6 = popt6_arr[:,1]/9
gamma6 = popt6_arr[:,2]/27

ind_Tc = np.argmin(abs(t_arr - 1))
gamma6_Tc = popt6_arr[ind_Tc,2]/27

ax4.plot(t_arr[:-1], -(alpha6)[:-1], label=r'$-\tilde\alpha(t)$')
# ax2.plot(t_arr, (t_arr - 1), 'k--')

line_beta = ax4.plot(t_arr, beta6/h.beta_const, label=r'$\tilde\beta(t)/\beta_0$')
# ax4.plot(t_arr, h.beta_B_norm(np.zeros_like(t_arr),0)/h.beta_const, ls='--', c=line_beta[0].get_c())
ax4.plot(t_arr, (5/3)/t_arr**2, ls='--')
# ax4.plot(t_arr, np.log(t_arr)/(t_arr-1), ls='--', c='b')

line_gamma6 = ax4.plot(t_arr, gamma6/gamma6_Tc, '.', label=r'$-\tilde\gamma(t)/\gamma(T_c)$')


# ax4.plot(t_arr, np.sqrt(-0.5*alpha6/beta6 ), 'k', label=r'$\sqrt{-\tilde{\alpha}/2\tilde{\beta}}$')
# ax4.plot(t_arr, bcs_gap_arr, 'b--', label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')


ax4.set_xlabel(r'$T/T_c$')

ax4.set_xlim(0,1.25)
ax4.set_ylim(0,20)

ax4.legend()
ax4.grid(True)

#%%

fig5, ax5 = plt.subplots()

def solve_quadratic(a, b, c):
    disc = np.sqrt(b**2 - 4*a*c)
    return (-b + disc)/(2*a), (-b - disc)/(2*a) 


# ax5.plot(t_arr, np.sqrt(-0.5*alpha6/beta6 )/np.sqrt(3), 'k', label=r'$\sqrt{-\tilde{\alpha}_6/2\tilde{\beta}_6}$')
# ax5.plot(t_arr, -0.5*alpha/beta/np.sqrt(3), 'r--', label=r'$\sqrt{-\tilde{\alpha}/2\tilde{\beta}}$')

# ax5.plot(t_arr, np.sqrt(3*alpha6*gamma6/beta6**2 ), 'g--', label=r'$3\alpha\gamma/\beta^2$')



approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6, 2*beta6, alpha6)

p = 22
delta_b_B = np.sum(h.delta_b_asarray(p, [1,2])) + np.sum(h.delta_b_asarray(p, [3,4,5]))/3 
delta_beta_B = h.beta_const * delta_b_B

approx_gap2_plus_sc, _  = solve_quadratic(3*gamma6, 2*(beta6 + t_arr*delta_beta_B), alpha6)

# approx_gap2_plus, approx_gap2_minus  = solve_quadratic(3*gamma6_Tc*np.ones_like(t_arr)/t_arr**3, 2*h.beta_A_norm(0*t_arr, 0)/t_arr**2, 
#                                                        t_arr - 1 - 0.5*(t_arr - 1)**2)


ax5.plot(t_arr, approx_gap2_plus/3, 
         label=r'$\Delta_+(\alpha,\beta,\gamma)/k_BT_c$')
ax5.plot(t_arr, approx_gap2_plus_sc/3, 
         label=r'$\Delta_+(\alpha,\beta^{\rm sc},\gamma)/k_BT_c$')
# ax5.plot(t_arr, np.sqrt(approx_gap2_minus)/np.sqrt(3), label=r'$\Delta_-$')


ax5.plot(t_arr, bcs_gap_arr**2, 'b--', label=r'$\Delta_{\rm BCS}(t)/k_BT_c$')


ax5.set_xlabel(r'$T/T_c$')
ax5.set_ylabel(r'$(\Delta/k_BT_c)^2$')

ax5.set_xlim(0,1.25)
ax5.set_ylim(0,5)

ax5.legend()
ax5.grid(True)