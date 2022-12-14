#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 18:02:02 2022

Stron coupling parameters comparator

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
import numpy.polynomial as nppoly


p_a = (0,34)
p_smooth = np.linspace(np.min(p_a), np.max(p_a), 200)

cols = ['r', 'g', 'b', 'c', 'm']
    
    
h.set_default('DEFAULT_SC_CORRS', 'WS15')
beta_WS15_arr = h.beta_norm_asarray(1, p_smooth)

for beta, col in zip(beta_WS15_arr, cols):
    plt.plot(p_smooth, beta, ls='--', c=col)


p_nodes = np.linspace(0,34,35)

h.set_default('DEFAULT_SC_CORRS', 'Choi')
beta_Choi_arr = h.beta_norm_asarray(1, p_nodes)

for beta, col in zip(beta_Choi_arr, cols):
    plt.plot(p_nodes, beta, 'x', c=col)

beta_Choi_data = h.get_interp_data_choi()

p_choi = beta_Choi_data[0]
beta_choi = beta_Choi_data[1]

my_poly_list = []
my_poly_order_list = [3, 8, 9, 7, 7]

for n, (beta, col, my_poly_order) in enumerate(zip(beta_choi, cols, my_poly_order_list)):
    plt.plot(p_choi, h.beta_const*(h.beta_norm_wc(n+1)+np.array(beta)), '.', c=col)
    my_poly_list.append(nppoly.Polynomial.fit(p_choi, beta_choi[n], my_poly_order))
    plt.plot(p_smooth, h.beta_const*(h.beta_norm_wc(n+1)+my_poly_list[n](p_smooth)), c=col)



plt.xlim(*p_a)
# plt.ylim(0.01, 0.03)

plt.legend()