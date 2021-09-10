#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:04:39 2021

@author: hindmars
"""

import numpy as np


# For method of linear interpolation
# From Regan, Wiman, Sauls arXiv:1908.04190 Table 1

p_nodes = range(0, 36, 2)

c1 = [-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275,
      -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, 
      -0.0402, -0.0413]

c2 = [-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, 
      -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, 
      -0.1583, -0.1645]

c3 = [-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, 
      -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, 
      -0.0267, -0.0268]

c4 = [-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, 
      -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, 
      -0.3388, -0.3518]

c5 = [-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, 
      -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, 
      -0.3717, -0.3815]

c_list = [c1, c2, c3, c4, c5]


# For polynomial fit method
# From Regan, Wiman, Sauls arXiv:1908.04190 Table 2
# Gives different values.  Typos, or pressure units?
# Check also Wiman thesis

# Only beta_1 at first.
a1 = [9.849e-3, -5.043e-2, 2.205e-2, -2.557e-2, 5.023e-2 -2.769e-2]



a_list = [a1]


# Functions

def delta_beta_norm_interp(p, n): 
    return np.interp(p, p_nodes, c_list[n-1])


def delta_beta_norm_poly(p, n): 
    
    if n==1:
        return np.poly1d(a_list[n-1], len(a_list[n-1]))(p)
    else:
        raise ValueError("only beta_1 implemented as yet")
        return


