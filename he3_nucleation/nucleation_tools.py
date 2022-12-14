#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 08:48:20 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt


def survival_prob(T, T_n, eta_n, T_0=0):
    t = (T - T_0)/(T_n - T_0)
    return np.exp(- t**(-eta_n))

def nuc_prob_dist(T, T_n, eta_n, norm, T_0=0):
    t = (T - T_0)/(T_n - T_0)
    return norm*eta_n * t**(- eta_n - 1) * np.exp(- t**(-eta_n))/(T_n - T_0)

def nuc_prob(T, T_n, eta_n, T_0=0):
    return 1 - survival_prob(T, T_n, eta_n, T_0)


