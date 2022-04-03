#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 19:09:45 2022

Calculates and plots the strong coupling corrections adjustment factor 
required to match the experimental T_AB.

Wiman thesis 2018 SC corrections, Figure 8.7 (scanned with WebPlotDigitizer)

The correction factor $f(p)$ is applied as
$$
\beta_a = \beta^{\rm wc} ( b_a + t e^f(p)  \Delta b_a^{\rm sc}(p) ),
$$
where $t$ is the reduced temperature, $p$ is the pressure, $b_a = (-1,2,2,2,-2)$, 
and $\Delta b_a^{\rm sc}(p)$ is the theoretical strong coupling correction.


@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
from scipy.optimize import curve_fit
# import numpy.polynomial as npp


h.DEFAULT_SC_CORRS = "Wiman-thesis"

p = np.linspace(h.p_pcp_bar,34,100)

# %%Fudge index
q = np.log(h.tAB_expt(p))/np.log(h.tAB(p))

plt.figure()
plt.plot(p, q, label='Index correction')

plt.xlabel(r'$p/bar$')
plt.ylabel(r'$q = \ln ( t^{\rm ex}_{AB})/\ln(t^{\rm th}_{AB})$')

plt.grid()
plt.legend()
plt.title("SC correction $(T/T_c)^{q(p)}$ - " + h.DEFAULT_SC_CORRS)


#%% Fudge factor

def fun(x, a0, a1, a2):
    return a0 + a1*x + a2*x**2

def fun2(x, a1, a2):
    return a1*x + a2*x**2

f = np.log(h.tAB_expt(p)/h.tAB(p))

plt.figure()

plt.plot(p, f, label='Correction exponent')

popt, copt = curve_fit(fun2, p[f>0] - h.p_pcp_bar, f[f>0], p0=(0.1,-0.1) )

plt.plot(p, fun2(p- h.p_pcp_bar, *popt) , 'r--',
          # label=r'$f \simeq {:3f} + {:.3f}\Delta p + ({:.3f})\Delta p^2$'.format(*popt))#, 
          label=r'$f \simeq {:.3g}\Delta p + ({:.3g})\Delta p^2$'.format(*popt))#, 


plt.xlabel(r'$p/bar$')
plt.ylabel(r'$f = \ln ( t^{\rm ex}_{AB}/t^{\rm th}_{AB})$')

plt.grid()
plt.legend()
plt.title("SC correction $(T/T_c) e^{f(p)}$ - " + h.DEFAULT_SC_CORRS)


plt.figure()


