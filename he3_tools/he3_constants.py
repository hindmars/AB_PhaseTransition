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

cphy = c.physical_constants

# Numerical constants
zeta3 = sp.zeta(3)
beta_const = 7 * zeta3/(80 * c.pi**2)
xiGL_const = np.sqrt(7 * zeta3 / 20)

# Physical constants
# Helium 3 mass in u
mhe3_u = 3.0160293 
mhe3_kg = mhe3_u * cphy["atomic mass constant"][0]
kB = c.k
R = c.R
N_A = c.N_A
hbar = c.hbar

a_bcs = 3.235 # Exponent for fit to BCS gap
delta_bcs0 = np.pi * np.exp(-np.euler_gamma)

beta_norm_wc_list = [-1, 2, 2, 2, -2]

