#!/usr/bin/env, python3
#, -*-, coding:, utf-8, -*-
"""
Created on Wed Aug 10 08:59:45 2022

BCS gap calculator

@author: hindmars
"""

import numpy as np
# import matplotlib.pyplot as plt
import scipy.integrate as scint
# import scipy.optimize as scopt
import scipy.special as scspe

from sympy.solvers import solve
from sympy import symbols

x = symbols('x')

# import he3_tools as h

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
bcs_const = np.sqrt(8*np.pi**2/(7*scspe.zeta(3)))

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

def free_energy(D_in, t):
    """
    D_in: normlaised eigenvalue(s) of 2x2 gap matrix Delta
    
    Single D: eigenvalues of the gap matrix assumed equal
    
    Pair of Ds, evaluate each and add

    Future
    - ? Give vector d and compute eigenvalues.
    
    """
    
    if isinstance(D_in, int) or isinstance(D_in, float):

        D = float(D_in)
    
        integral_1 = 2*scint.quad(full_integrand_1, xmin, xmax, args=(D, t))[0]
        integral_2 = 2*scint.quad(full_integrand_2, xmin, xmax, args=(D, t))[0]
        
        factor = 2.0
        
    elif len(D_in) == 2:

        integral_1 = 0.0
        integral_2 = 0.0

        for D in D_in:
            integral_1 += 2*scint.quad(full_integrand_1, xmin, xmax, args=(D, t))[0]
            integral_2 += 2*scint.quad(full_integrand_2, xmin, xmax, args=(D, t))[0]
        
        factor = 1.0
    
    return 1.5*factor*((integral_1 + np.log(t))*D**2 - (integral_2/t + 0.5*D**2) )

def gap_eqn(D, t):
    def fun(x):
        # return - integrand_1(x, 0, 1) + integrand_1(x, D, t) 
        return 2*full_integrand_1(x, D, t) 
    res, _ = scint.quad(fun, xmin, xmax)  
    return res + np.log(t)

def gap_GL(t): 
    return bcs_const * np.sqrt(1-t)

def gap_eqn_A_phase(D, t):

    def gap_eqn_mu(mu):
        return D**2 * (1 - mu**2) *gap_eqn(D * np.sqrt(1 - mu**2), t)
    
    res, _ = scint.quad(gap_eqn_mu, -1, 1)

    return 0.5*res

def free_energy_A_phase(D, t):
    
    def free_energy_mu(mu):
        
        return free_energy(D*np.sqrt(1 - mu**2), t)
    
    res, _ = scint.quad(free_energy_mu, -1, 1)

    return 0.5*res

def quartic_symmetric(x, *a):
    return a[0]*x**2 + a[1]*x**4

def sextic_symmetric(x, *a):
    return a[0]*x**2 + a[1]*x**4 + a[2]*x**6

def poly_in_squares(x, *a):
    
    y = x * 0.0
    
    for n, aa in enumerate(a):
        y += aa*x**(2*n + 2)
        
    return y

def solve_quadratic(a,b,c):
    disc = np.sqrt(b**2 - 4*a*c)
    return 0.5*(-b + disc)/a, 0.5*(-b - disc)/a

def solve_cubic_arr(*a_arrs):
    
    sols = []
    for a in zip(*a_arrs):
        sols.append(solve(a[0] + a[1]*x + a[2]*x**2 + a[3]*x**3, x))
        
    return np.array(sols)

#%% From https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 5))),
     ('densely dotted',        (0, (1, 1))),

     ('long dash', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

linestyle_dict = {}

for name, spec in linestyle_tuple:
    linestyle_dict[name] = spec

ls_dict = linestyle_dict

color_list = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#999999']

