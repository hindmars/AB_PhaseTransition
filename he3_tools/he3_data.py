#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 15:01:24 2022

Contains experimental data and functions to return values of material parameters 
and other important quantities.

@author: hindmars
"""

import numpy as np
import numpy.polynomial as nppoly
import he3_constants as h3c

def convert_b_to_db(b_list, n):
    db_list = []
    for b in b_list:
        db_list.append(b - h3c.b_wc_list[n-1])
    return db_list


# Construct dictionaries of model strong coupling corrections.
# From Regan, Wiman, Sauls arXiv:1908.04190 Table 1

def get_interp_data_rws19():
    p_nodes_beta = range(0, 36, 2)
    
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

    return [p_nodes_beta, c_list]

def get_interp_data_wiman_thesis():
    """From Figure 8, grabbed using Web Plot Digitiser.
    """
    
    p_nodes = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 
               22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0]
    c1 = [-0.004830918, -0.008756039, -0.009057971, -0.011171498, -0.013586957, 
          -0.016304348, -0.01781401, -0.019021739, -0.021437198, -0.023852657, 
          -0.025664251, -0.027475845, -0.028683575, -0.030495169, -0.031702899, 
          -0.033514493, -0.032689211, -0.033977456]
    c2 = [-0.025362319, -0.039855072, -0.054951691, -0.068236715, -0.080012077, 
          -0.090881643, -0.101751208, -0.111714976, -0.120471014, -0.128623188, 
          -0.136775362, -0.144625604, -0.151268116, -0.157306763, -0.162741546, 
          -0.168176329, -0.173007246, -0.177536232]
    c3 = [-0.014492754, -0.020833333, -0.023852657, -0.026871981, -0.02928744, 
          -0.032306763, -0.035929952, -0.03955314, -0.041364734, -0.041364734, 
          -0.040458937, -0.041364734, -0.041666667, -0.041968599, -0.041968599, 
          -0.042572464, -0.041968599, -0.041062802]
    c4 = [-0.026570048, -0.036533816, -0.044987923, -0.054649758, -0.065519324, 
          -0.077294686, -0.089070048, -0.102355072, -0.115640097, -0.127717391, 
          -0.139794686, -0.155495169, -0.172101449, -0.189915459, -0.210144928, 
          -0.233091787, -0.257850242, -0.282608696]
    c5 = [-0.079710145, -0.129307568, -0.174214976, -0.213768116, -0.248792271, 
          -0.282608696, -0.310386473, -0.33544686, -0.358997585, -0.380434783, 
          -0.400664251, -0.415157005, -0.4272343, -0.4375, -0.443538647, 
          -0.446557971, -0.444746377, -0.443236715]

    c_list = [c1, c2, c3, c4, c5]

    return [p_nodes, c_list]

def get_interp_data_choi():
    p_choi = range(0, 35, 1)
    b1_choi = [-0.97, -0.97, -0.97, -0.98, -0.98, -0.98, -0.98, -0.98, -0.98, -0.99, -0.99, -0.99, -0.99, -0.99,
               -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.01, -1.01, -1.01, -1.01, -1.01, -1.01, -1.02, -1.02,
               -1.02,  -1.02, -1.02, -1.03, -1.03, -1.03, -1.03]
    b2_choi = [1.89, 1.94, 1.96, 1.99, 1.99, 1.99, 1.99, 1.98, 1.98, 1.98, 1.97, 1.97, 1.96, 1.95, 1.95, 1.95, 1.95,
               1.94,1.94, 1.93, 1.94, 1.94, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93,
               1.93]
    b3_choi = [2.10, 1.96, 1.86, 1.81, 1.76, 1.74, 1.72, 1.70, 1.70, 1.69, 1.69, 1.70, 1.69, 1.69, 1.70, 1.72, 1.73, 1.72,
               1.73, 1.72, 1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 1.73, 1.74, 1.73, 1.73, 1.72, 1.73, 1.73, 1.73,
               1.73]
    b4_choi = [1.85, 1.72, 1.63, 1.56, 1.52, 1.48, 1.46, 1.44, 1.42, 1.41, 1.40, 1.39, 1.39, 1.39, 1.38, 1.35, 1.34,
               1.33,1.32, 1.33, 1.31, 1.29, 1.29, 1.29, 1.28, 1.28, 1.27, 1.26, 1.26, 1.26, 1.26, 1.25, 1.25, 1.25,
               1.25]
    b5_choi = [-1.84, -1.82, -1.81, -1.81, -1.81, -1.81, -1.82, -1.82, -1.83, -1.84, -1.85, -1.86, -1.87, -1.88, -1.89, -1.89, -1.90,
               -1.90, -1.91, -1.92, -1.93, -1.93, -1.94, -1.95, -1.96, -1.97, -1.97, -1.98, -1.99, -2.00, -2.01, -2.02, -2.02, -2.03,
               -2.03]
    
    c1 = convert_b_to_db(b1_choi, 1)
    c2 = convert_b_to_db(b2_choi, 2)
    c3 = convert_b_to_db(b3_choi, 3)
    c4 = convert_b_to_db(b4_choi, 4)
    c5 = convert_b_to_db(b5_choi, 5)
    
    c_choi_list = [c1, c2, c3, c4, c5]

    return [p_choi, c_choi_list]

def get_poly_ws15():
    # polynomial coefficients from J. J. Wiman ‚ J. A. Sauls prb 92, 144515 (2015)
    # \Delta{\beta_i^{sc}}/|beta_1^{wc}| = a_n^{i} p^{n}, p in bar
    a1sc = [3.070e-2, -2.081e-3, 2.133e-5, -4.189e-7]
    a2sc = [-1.074e-1, 5.412e-2, -1.081e-2, 1.025e-3, -5.526e-5, 1.722e-6, -2.876e-8, 1.991e-10]
    a3sc = [1.038e-1, -1.752e-1, 3.488e-2, -4.243e-3, 3.316e-4, -1.623e-5, 4.755e-7, -7.587e-9, 5.063e-11]
    a4sc = [-1.593e-1, -1.350e-1, 1.815e-2, -1.339e-3, 5.316e-5, -1.073e-6, 8.636e-9]
    a5sc = [1.610e-1, 2.263e-2, -4.921e-3, 3.810e-4, -1.529e-5, 3.071e-7, -2.438e-9]
    return [nppoly.Polynomial(a1sc), nppoly.Polynomial(a2sc), nppoly.Polynomial(a3sc), 
            nppoly.Polynomial(a4sc), nppoly.Polynomial(a5sc)]

def get_poly_choi_this():
    # Polynomial fits to Choi et al data, same orders as Wiman-Sauls 2015
    beta_Choi_data = get_interp_data_choi()

    p_choi = beta_Choi_data[0]
    beta_choi = beta_Choi_data[1]

    my_poly_list = []
    my_poly_order_list = [3, 8, 9, 7, 7]

    for n, (beta, my_poly_order) in enumerate(zip(beta_choi, my_poly_order_list)):
        my_poly_list.append(nppoly.Polynomial.fit(p_choi, beta_choi[n], my_poly_order))

    return my_poly_list

def get_poly_rws19_this():
    """Polynomial fits to Regan-Wiman-Sauls 2019 quoted Wiman beta data, same orders 
    as they give, but fitted directly."""
    beta_rws19_data = get_interp_data_rws19()

    p_rws19 = beta_rws19_data[0]
    beta_rws19 = beta_rws19_data[1]

    my_poly_list = []
    my_poly_order_list = [5, 5, 5, 5, 5]

    for n, (beta, my_poly_order) in enumerate(zip(beta_rws19, my_poly_order_list)):
        my_poly_list.append(nppoly.Polynomial.fit(p_rws19, beta_rws19[n], my_poly_order))

    return my_poly_list


dbeta_data_dict = { "RWS19" : get_interp_data_rws19(),
                   "RWS19-interp" : get_interp_data_rws19(),
                   "RWS19-poly" : get_poly_rws19_this(),
                "Wiman-thesis" : get_interp_data_wiman_thesis(),
                "Choi-interp" : get_interp_data_choi(),
                "Choi-poly" : get_poly_choi_this(),
                "WS15" : get_poly_ws15(),
                "WS15-poly" : get_poly_ws15()
                }

######################################################################################
### Interpolation data for other material parameters, from Greywall 1986 via RWS19 ###
######################################################################################
##'''
#
p_nodes = list(range(0, 36, 2))

Tc_data_mK = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 
      2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]

# particle_density, nm^{-3}
np_data = [16.28, 17.41, 18.21, 18.85, 19.34, 19.75, 20.16, 20.60, 21.01, 21.44, 21.79, 
       22.96, 22.36, 22.54, 22.71, 22.90, 23.22, 23.87]

# Effective mass ratio
mstar_m_data = [2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 
           5.02, 5.18, 5.34, 5.50, 5.66, 5.8]

# Fermi velocity, m/s
vf_data = [59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 
      37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23]

# Coherence length at T=0, nm
xi0_data = [77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 
       19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76]

# Landau parameter $F_0^a$
F0a_data = [-0.7226, -0.7317, -0.7392, -0.7453, -0.7503, -0.7544, -0.7580, -0.7610, 
            -0.7637, -0.7661, -0.7684, -0.7705, -0.7725, -0.7743, -0.7758, -0.7769, 
            -0.7775, -0.7775]
F0a_poly = nppoly.Polynomial.fit(p_nodes, F0a_data, 5)

# For polynomial fit method Tc
# From Regan, Wiman, Sauls arXiv:1908.04190 Table 2
# Doesn't work at the moment, some unit miunserstanding or misprint?
# a1 = [-9.849e-3, -5.043e-2, 2.205e-2, -2.557e-2, 5.023e-2 -2.769e-2]
# a_list = [a1[::-1]] # polyfit wants highest power first
# b1_poly = np.polynomial.Polynomial(a1)

###############################################################################################################
##########            Greywall, PLTS temperature scales polynomial coefficients                 ###############
###############################################################################################################


# From Greywall 1986
T_pcp_mK = 2.273
p_pcp_bar = 21.22
p_A_bar = 34.338
# Greywall 1986 Tc polynomial fit cofficients and constructor
aTc_G = [0.92938375, 0.13867188, -0.69302185e-2, 0.25685169e-3, -0.57248644e-5, 0.5301091e-7]
Tc_poly_Greywall = nppoly.Polynomial(aTc_G)
# Greywall 1986 TAB polynomial fit cofficients and constructor
aTAB_G = [T_pcp_mK, -0.10322623e-1, -0.53633181e-2, 0.83437032e-3, -0.61709783e-4,  0.17038992e-5]
TAB_poly_Greywall = nppoly.Polynomial(aTAB_G)

# Melting curve poly coeffs (Greywall) 0.9 - 5.6 mK. 
aPmelt = [-0.19632970e-1, 0.61880268e-1, -0.78803055e-1, 0.13050600, -0.43519381e-1, 
          0.13752791e-3, -0.17180436e-6, -0.220939060e-9, 0.85450245e-12] 
Pmelt_poly_Greywall = nppoly.Polynomial(aPmelt)

# Greywall scale to PLTS scale, 6 order polynomial coefficients (Parpia et al 2022)
# Temperature 0.9 mK to 5.6 mK
GtoPLTS6 = [-0.14265343150487, 1.2810635032153, -0.22689947807354, 0.084337673002034, 
            -0.016928990685839, 0.0017611612884063, -7.4461876859237e-5]
GtoPLTS_low_poly = nppoly.Polynomial(GtoPLTS6)

# Greywall scale to PLTS scale, 9 order polynomial coefficients (Parpia et al 2022)
# Temperature 5.6 mK to 100 mK
GtoPLTS9 = [0.020353327019475, 0.96670033496024, 0.0019559314169033, -9.5551084662924e-5, 
            3.2167457655106e-6, -7.0097586342143e-8, 9.6909878738352e-10,
            -8.2126513949290e-12, 3.8886762300964e-14, -7.8713540127550e-17]
GtoPLTS_high_poly = nppoly.Polynomial(GtoPLTS9)

# Invert the transformation
def invert_poly(poly, n, x_range):
    x = np.linspace(x_range[0], x_range[-1], 100)
    return nppoly.Polynomial.fit(poly(x), x, n)
    
PLTStoG9_low_poly = invert_poly(GtoPLTS_low_poly, 9, [0.9, 5.6])
PLTStoG9_high_poly = invert_poly(GtoPLTS_high_poly, 9, [5.6, 100])

# PLTS scale coefficients for polynomial method for Tc and T_AB, from Tian, Smith, Parpia et al 2022
d_c = np.array([0.90972399274531, 0.14037182852625, -0.0074017331747577, 2.8617547367067e-4,-6.5064429600510e-6, 6.0754459040296e-8])
c_AB = np.array([-26.864685876026, 5.2647866128370, -0.37617826876151, 0.013325635880953, -2.3510107585468e-4, 1.6519539175010e-6])

Tc_poly_PLTS = nppoly.Polynomial(d_c)
TAB_poly_PLTS = nppoly.Polynomial(c_AB)

# PLTS temperature along melting curve, 5.6 - 100 mK, p - pA in millibar
b_melt_hi = np.array([
    2.5257068036689, 
    -0.024555787070591, 
    -1.7994349872600e-6, 
    -6.0072773787192e-9, 
    -5.8885557756054e-12, 
    -3.9102041149206e-15, 
    -1.6605359484626e-18, 
    -4.3792852104458e-22, 
    -6.5042399327047e-26, 
    -4.1677548758514e-30] )

Tmelt_poly_PLTS_hi = nppoly.Polynomial(b_melt_hi) 

# PLTS temperature along melting curve, 0.9 - 5.6 mK, p - pA in millibar
b_melt_lo = np.array([
    2.4442188707375, 
    -0.026522783446809, 
    -2.7665556176467e-5, 
    -2.2800036357244e-7, 
     3.9890194953355e-10, 
     1.2845430171276e-11, 
    -6.9521369379387e-13, 
    -1.5658128424388e-14, 
    -1.1687750824147e-16, 
    -3.0194721850282e-19] )

Tmelt_poly_PLTS_lo = nppoly.Polynomial(b_melt_lo) # convert to bar


