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

r""" The BCS equation for the gap $\Delta$ is (Vollhardt and Woelfle 3.46e)
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

xmax = 1e2 # Cut-off on integrations
xmin = 1e-2
const = np.sqrt(8*np.pi**2/(7*scspe.zeta(3)))

def integrand(x, D, t=1):
    e = np.sqrt(x**2 + D**2)
    return np.tanh(0.5*e/t)/e

def full_integrand(x, D, t):
    return integrand(x, D, t) - integrand(x, 0) 

def gap_eqn(D, t):
    def fun(x):
        return full_integrand(x, D, t)
    res, _ = scint.quad(fun, xmin, xmax)
    return res

def gap_GL(t): 
    return const * np.sqrt(1-t)

t_arr = np.linspace(0.002, 1, 500, endpoint=False)

D_arr = np.array([float(scopt.fsolve(gap_eqn, 0.1, args=(t,))) for t in t_arr])

#%%
# Fit the gap

def gap_fit1(t,a):
    """
    
    Best fit a [3.23485079]
    
    """
    return (1 - t**a)**0.5

def gap_fit2(t,a,b):
    """
    
    Best fit [a,b] [3.48540465, 0.541188  ]
    
    """
    return (1 - t**a)**b

def gap_fit3(t, *pars):
    """

    Best fit pars: [-2.06196056e-03,  3.13789240e+00, -3.11236328e+00,  9.61261152e-01]    

    """
    x = 1-t
    s = 0
    for n, p in enumerate(pars):
        s += p*x**(n+1)
    return s**0.5

def gap_fit4(t, *pars):
    """

    Best fit pars: 

    """
    # x = 1-t
    s = 1 - t**2
    # for n, p in enumerate(pars[:-2]):
    #     s += p*x**(n+2)
    s = s/(pars[0] + pars[1]*t)
    return np.abs(s)**0.5

def make_lab3(*pars):
    lab3 = ''
    for n, p in enumerate(pars):
        if n== 0:
            lab3 += r'$[{:.3g}x '.format(pars[n])
        else:
            lab3 += '+ ({:.3g})x^{:d} '.format(pars[n], n+1)
        
    lab3 += ']^{1/2}$'
    return lab3
        
def make_lab4(*pars):
    lab4 = ''
    # lab4 = r'$[x '
    # for n, p in enumerate(pars[:-2]):
    #     # if n== 0:
    #     #     lab4 += r'$[{:.3g} '.format(pars[n])
    #     # else:
    #     lab4 += '+ ({:.3g})x^{:d} '.format(pars[n], n+2)
        
    # lab4 += r']^{{1/2}}/[({:.3g}) + ({:.3g})t]^{{1/2}}$'.format(*pars[-2:])

    lab4 += r'$[ 1 - t^2]^{{1/2}}/[({:.3g}) + ({:.3g})t]^{{1/2}}$'.format(*pars[-2:])
    return lab4

#%%
popt1, pcov1 = scopt.curve_fit(gap_fit1, t_arr, D_arr/D_arr[0], p0=(4,))
popt2, pcov2 = scopt.curve_fit(gap_fit2, t_arr, D_arr/D_arr[0], p0=(4,0.5))
popt3, pcov3 = scopt.curve_fit(gap_fit3, t_arr, D_arr/D_arr[0], p0=(1,1,1))
popt4, pcov4 = scopt.curve_fit(gap_fit4, t_arr, D_arr/D_arr[0], p0=(1,-1))

print(popt1, popt2, popt3,popt4)

#%%
# From Muehlschleger 1959
muehl_data = np.array([1.00, 0.0000, 0.98, 0.2436, 0.96, 0.3416, 0.94, 0.4148, \
                       0.92, 0.4749, 0.90, 0.5263, 0.88, 0.5715, 0.86, 0.6117, \
                       0.84, 0.6480, 0.82, 0.6810, 0.80, 0.7110, 0.78, \
                       0.7386, 0.76, 0.7640, 0.74, 0.7874, 0.72, 0.8089, 0.70, \
                       0.8288, 0.68, 0.8474, 0.66, 0.8640, 0.64, 0.8796, 0.62, \
                       0.8939, 0.60, 0.9070, 0.58, 0.9190, 0.56, 0.9299, 0.54, \
                       0.9399, 0.52, 0.9488, 0.50, 0.9569, 0.48, 0.9641, 0.46, \
                       0.9704, 0.44, 0.9760, 0.42, 0.9809, 0.40, 0.9850, 0.38, \
                       0.9885, 0.36, 0.9915, 0.34, 0.9938, 0.32, 0.9957, 0.30, \
                       0.9971, 0.28, 0.9982, 0.26, 0.9989, 0.24, 0.9994, 0.22, \
                       0.9997, 0.20, 0.9999, 0.18, 1.0000, 0.16, 1.0000, 0.14, 
                       1.0000])

muehl_t = muehl_data[0::2]
muehl_delta = muehl_data[1::2]
#%%

plt.plot(t_arr, gap_GL(t_arr), 'k--', label='GL gap')
plt.plot(t_arr, D_arr, label='BCS gap')
# plt.plot(t_arr, D_arr[0]*gap_fit1(t_arr, *popt1), 
#          label=r'$\Delta_0(1 - t^{{ {:.3f} }})^{{ 0.5 }}$'.format(*popt1))
# plt.plot(t_arr, D_arr[0]*gap_fit2(t_arr, *popt2), 
#           label=r'$\Delta_0(1 - t^{{ {:.3f} }})^{{ {:.3f} }}$'.format(*popt2))
# plt.plot(t_arr, D_arr[0]*gap_fit3(t_arr, *popt3), 
#           label=make_lab3(*popt3))
plt.plot(t_arr, D_arr[0]*gap_fit4(t_arr, *popt4), 
          label=make_lab4(*popt4))
plt.plot(muehl_t, D_arr[0]*muehl_delta,'k.', label=r'Muehlschlegel 1959')

plt.xlabel(r'$T/T_c$')
plt.ylabel(r'$\Delta_T/T_c$')
plt.xlim(0,1)
plt.ylim(0,3)
plt.legend()
plt.grid(True)

plt.savefig('bcs_gap_and_fit.pdf')


