#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 15:16:18 2022

Useful constants for he3_tools.

@author: hindmars
"""
import scipy.special as sp
import scipy.constants as c
import numpy as np

ibm_colorblind_list = ['#648fff', '#785ef0', '#dc267f', '#fe6100', '#ffb000', '#000000', '#ffffff']
alt_colorblind_list = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#333333']
wong_colourblind_list = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']


cphy = c.physical_constants

# Numerical constants
zeta3 = sp.zeta(3)
beta_const = 7 * zeta3/(80 * c.pi**2)
xiGL_const = np.sqrt(7 * zeta3 / 20)
b_wc_list = [-1, 2, 2, 2, -2]

# Physical constants
kB = c.k
R = c.R
N_A = c.N_A
hbar = c.hbar

# Helium 3 constants

mhe3_u = 3.0160293 
mhe3_kg = 3.0160293 * cphy['atomic mass constant'][0]
mu0he3_J_T = cphy["helion mag. mom."][0]

# Zeeman (quadratic) energy constant, evaaluated at 1 T and 1 mK
gH0 = 20 * beta_const * (mu0he3_J_T/kB)**2 * (1/1e-3)**2
# PHA asymmetry parameter for linear magnetic energy constant (see Sauls and Sharma 2003)
lambda_A1 = 60.1e-3 # mK/T

# BCS constants
a_bcs = 3.235 # Exponent for fit to BCS gap
delta_bcs0 = np.pi * np.exp(-np.euler_gamma)


