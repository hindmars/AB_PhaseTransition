#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:04:39 2021

@author: hindmars;  
"""

import numpy as np
# import scipy.linalg as sl
import pandas as pd
import numpy.polynomial as nppoly

import he3_bases as h3b
import he3_constants as h3c
import he3_data as h3d
# import he3_matrix as h3m

# import scipy.special as sp
# import scipy.constants as c

# cphy = c.physical_constants

# # Numerical constants
# zeta3 = sp.zeta(3)
# beta_const = 7 * zeta3/(80 * np.pi**2)
# xiGL_const = np.sqrt(7 * zeta3 / 20)

# # Physical constants
# # Helium 3 mass in u
# mhe3_u = 3.0160293 
# mhe3_kg = mhe3_u * cphy["atomic mass constant"][0]
# kB = c.k
# R = c.R
# N_A = c.N_A
# hbar = c.hbar

# a_bcs = 3.235 # Exponent for fit to BCS gap
# delta_bcs0 = np.pi * np.exp(-np.euler_gamma)

# beta_norm_wc_list = [-1, 2, 2, 2, -2]


SET_T_SCALE= {"Greywall", "PLTS"}
# DEFAULT_T_SCALE="Greywall" 
DEFAULT_T_SCALE="PLTS" 

sc_corrs_interp = {"RWS19", "RWS19-interp", "Wiman-thesis", "Choi-interp"}
sc_corrs_poly = {"RWS19-poly", "Choi-poly", "WS15", "WS15-poly"}
# SET_SC_CORRS= {"RWS19", "Wiman-thesis", "Choi-interp", "WS15", "Choi-poly"}
SET_SC_CORRS= sc_corrs_interp.union(sc_corrs_poly)

DEFAULT_SC_CORRS="RWS19"
# DEFAULT_SC_CORRS="Wiman-thesis"
# DEFAULT_SC_CORRS="Choi"

# Do we want to adjust the SC corrections to get TAB right?
SET_SC_ADJUST = {True, False}
DEFAULT_SC_ADJUST=False
SET_SC_ADJUST_MIN_P={0, "p_pcp_bar"}
DEFAULT_SC_ADJUST_MIN_P="p_pcp_bar" # Should set this as p_pcp_bar
# Polynomial for adjusting SC corrections to fit TAB data
# Gets redefined later on if DEFAULT_SC_ADJUST==True.
sc_corr_adj_pol = nppoly.Polynomial([0])

SET_ALPHA_TYPE = {"GL", "BCS"}
DEFAULT_ALPHA_TYPE = "GL"

def report_setting(name):
    xval = globals()[name]
    print("he3_tools:", type(xval), "variable " + name + " set to", xval)

def set_default(name, xval):
    set_name = name.replace("DEFAULT_", "SET_")
    if xval in globals()[set_name]:
        globals()[name] = xval
    else:
        raise ValueError("error: " + xval + " not in " + set_name)
    report_setting(name)
    return 

def beta_norm_wc(n):
    """
    Parameters
    ----------
    n : int
        Label for GL bulk free energy beta parameter.

    Returns
    -------
    b : float (or python default, usually is double)
        Normalised beta, in weak coupling app.

    """
    b = np.nan
    if n in [1,2,3,4,5]:
        b = h3c.beta_norm_wc_list[n-1]
    else:
        raise ValueError("beta_norm_wc: n should be 1, 2, 3, 4 , or 5")
    return b

# def convert_b_to_db(b_list, n):
#     db_list = []
#     for b in b_list:
#         db_list.append(b - beta_norm_wc(n))
#     return db_list

# # Construct dictionaries of model strong coupling corrections.
# # From Regan, Wiman, Sauls arXiv:1908.04190 Table 1

# def get_interp_data_rws19():
#     p_nodes_beta = range(0, 36, 2)
    
#     c1 = [-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275,
#           -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, 
#           -0.0402, -0.0413]
    
#     c2 = [-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, 
#           -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, 
#           -0.1583, -0.1645]
    
#     c3 = [-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, 
#           -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, 
#           -0.0267, -0.0268]
    
#     c4 = [-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, 
#           -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, 
#           -0.3388, -0.3518]
    
#     c5 = [-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, 
#           -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, 
#           -0.3717, -0.3815]

#     c_list = [c1, c2, c3, c4, c5]

#     return [p_nodes_beta, c_list]

# def get_interp_data_wiman_thesis():
#     """From Figure 8, grabbed using Web Plot Digitiser.
#     """
    
#     p_nodes = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 
#                22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0]
#     c1 = [-0.004830918, -0.008756039, -0.009057971, -0.011171498, -0.013586957, 
#           -0.016304348, -0.01781401, -0.019021739, -0.021437198, -0.023852657, 
#           -0.025664251, -0.027475845, -0.028683575, -0.030495169, -0.031702899, 
#           -0.033514493, -0.032689211, -0.033977456]
#     c2 = [-0.025362319, -0.039855072, -0.054951691, -0.068236715, -0.080012077, 
#           -0.090881643, -0.101751208, -0.111714976, -0.120471014, -0.128623188, 
#           -0.136775362, -0.144625604, -0.151268116, -0.157306763, -0.162741546, 
#           -0.168176329, -0.173007246, -0.177536232]
#     c3 = [-0.014492754, -0.020833333, -0.023852657, -0.026871981, -0.02928744, 
#           -0.032306763, -0.035929952, -0.03955314, -0.041364734, -0.041364734, 
#           -0.040458937, -0.041364734, -0.041666667, -0.041968599, -0.041968599, 
#           -0.042572464, -0.041968599, -0.041062802]
#     c4 = [-0.026570048, -0.036533816, -0.044987923, -0.054649758, -0.065519324, 
#           -0.077294686, -0.089070048, -0.102355072, -0.115640097, -0.127717391, 
#           -0.139794686, -0.155495169, -0.172101449, -0.189915459, -0.210144928, 
#           -0.233091787, -0.257850242, -0.282608696]
#     c5 = [-0.079710145, -0.129307568, -0.174214976, -0.213768116, -0.248792271, 
#           -0.282608696, -0.310386473, -0.33544686, -0.358997585, -0.380434783, 
#           -0.400664251, -0.415157005, -0.4272343, -0.4375, -0.443538647, 
#           -0.446557971, -0.444746377, -0.443236715]

#     c_list = [c1, c2, c3, c4, c5]

#     return [p_nodes, c_list]

# def get_interp_data_choi():
#     p_choi = range(0, 35, 1)
#     b1_choi = [-0.97, -0.97, -0.97, -0.98, -0.98, -0.98, -0.98, -0.98, -0.98, -0.99, -0.99, -0.99, -0.99, -0.99,
#                -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.01, -1.01, -1.01, -1.01, -1.01, -1.01, -1.02, -1.02,
#                -1.02,  -1.02, -1.02, -1.03, -1.03, -1.03, -1.03]
#     b2_choi = [1.89, 1.94, 1.96, 1.99, 1.99, 1.99, 1.99, 1.98, 1.98, 1.98, 1.97, 1.97, 1.96, 1.95, 1.95, 1.95, 1.95,
#                1.94,1.94, 1.93, 1.94, 1.94, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93, 1.93,
#                1.93]
#     b3_choi = [2.10, 1.96, 1.86, 1.81, 1.76, 1.74, 1.72, 1.70, 1.70, 1.69, 1.69, 1.70, 1.69, 1.69, 1.70, 1.72, 1.73, 1.72,
#                1.73, 1.72, 1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 1.73, 1.74, 1.73, 1.73, 1.72, 1.73, 1.73, 1.73,
#                1.73]
#     b4_choi = [1.85, 1.72, 1.63, 1.56, 1.52, 1.48, 1.46, 1.44, 1.42, 1.41, 1.40, 1.39, 1.39, 1.39, 1.38, 1.35, 1.34,
#                1.33,1.32, 1.33, 1.31, 1.29, 1.29, 1.29, 1.28, 1.28, 1.27, 1.26, 1.26, 1.26, 1.26, 1.25, 1.25, 1.25,
#                1.25]
#     b5_choi = [-1.84, -1.82, -1.81, -1.81, -1.81, -1.81, -1.82, -1.82, -1.83, -1.84, -1.85, -1.86, -1.87, -1.88, -1.89, -1.89, -1.90,
#                -1.90, -1.91, -1.92, -1.93, -1.93, -1.94, -1.95, -1.96, -1.97, -1.97, -1.98, -1.99, -2.00, -2.01, -2.02, -2.02, -2.03,
#                -2.03]
    
#     c1 = convert_b_to_db(b1_choi, 1)
#     c2 = convert_b_to_db(b2_choi, 2)
#     c3 = convert_b_to_db(b3_choi, 3)
#     c4 = convert_b_to_db(b4_choi, 4)
#     c5 = convert_b_to_db(b5_choi, 5)
    
#     c_choi_list = [c1, c2, c3, c4, c5]

#     return [p_choi, c_choi_list]

# def get_poly_ws15():
#     # polynomial coefficients from J. J. Wiman ‚ J. A. Sauls prb 92, 144515 (2015)
#     # \Delta{\beta_i^{sc}}/|beta_1^{wc}| = a_n^{i} p^{n}, p in bar
#     a1sc = [3.070e-2, -2.081e-3, 2.133e-5, -4.189e-7]
#     a2sc = [-1.074e-1, 5.412e-2, -1.081e-2, 1.025e-3, -5.526e-5, 1.722e-6, -2.876e-8, 1.991e-10]
#     a3sc = [1.038e-1, -1.752e-1, 3.488e-2, -4.243e-3, 3.316e-4, -1.623e-5, 4.755e-7, -7.587e-9, 5.063e-11]
#     a4sc = [-1.593e-1, -1.350e-1, 1.815e-2, -1.339e-3, 5.316e-5, -1.073e-6, 8.636e-9]
#     a5sc = [1.610e-1, 2.263e-2, -4.921e-3, 3.810e-4, -1.529e-5, 3.071e-7, -2.438e-9]
#     return [nppoly.Polynomial(a1sc), nppoly.Polynomial(a2sc), nppoly.Polynomial(a3sc), 
#             nppoly.Polynomial(a4sc), nppoly.Polynomial(a5sc)]

# def get_poly_choi_this():
#     # Polynomial fits to Choi et al data, same orders as Wiman-Sauls 2015
#     beta_Choi_data = get_interp_data_choi()

#     p_choi = beta_Choi_data[0]
#     beta_choi = beta_Choi_data[1]

#     my_poly_list = []
#     my_poly_order_list = [3, 8, 9, 7, 7]

#     for n, (beta, my_poly_order) in enumerate(zip(beta_choi, my_poly_order_list)):
#         my_poly_list.append(nppoly.Polynomial.fit(p_choi, beta_choi[n], my_poly_order))

#     return my_poly_list

# def get_poly_rws19_this():
#     """Polynomial fits to Regan-Wiman-Sauls 2019 quoted Wiman beta data, same orders 
#     as they give, but fitted directly."""
#     beta_rws19_data = get_interp_data_rws19()

#     p_rws19 = beta_rws19_data[0]
#     beta_rws19 = beta_rws19_data[1]

#     my_poly_list = []
#     my_poly_order_list = [5, 5, 5, 5, 5]

#     for n, (beta, my_poly_order) in enumerate(zip(beta_rws19, my_poly_order_list)):
#         my_poly_list.append(nppoly.Polynomial.fit(p_rws19, beta_rws19[n], my_poly_order))

#     return my_poly_list


# dbeta_data_dict = { "RWS19" : get_interp_data_rws19(),
#                    "RWS19-interp" : get_interp_data_rws19(),
#                    "RWS19-poly" : get_poly_rws19_this(),
#                 "Wiman-thesis" : get_interp_data_wiman_thesis(),
#                 "Choi-interp" : get_interp_data_choi(),
#                 "Choi-poly" : get_poly_choi_this(),
#                 "WS15" : get_poly_ws15(),
#                 "WS15-poly" : get_poly_ws15()
#                 }

# ######################################################################################
# ### Interpolation data for other material parameters, from Greywall 1986 via RWS19 ###
# ######################################################################################
# ##'''
# #
# p_nodes = list(range(0, 36, 2))

# Tc_data_mK = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 
#       2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]

# # particle_density, nm^{-3}
# np_data = [16.28, 17.41, 18.21, 18.85, 19.34, 19.75, 20.16, 20.60, 21.01, 21.44, 21.79, 
#        22.96, 22.36, 22.54, 22.71, 22.90, 23.22, 23.87]

# # Effective mass ratio
# mstar_m_data = [2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 
#            5.02, 5.18, 5.34, 5.50, 5.66, 5.8]

# # Fermi velocity, m/s
# vf_data = [59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 
#       37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23]

# # Coherence length at T=0, nm
# xi0_data = [77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 
#        19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76]

# # Landau parameter $F_0^a$
# F0a_data = [-0.7226, -0.7317, -0.7392, -0.7453, -0.7503, -0.7544, -0.7580, -0.7610, 
#             -0.7637, -0.7661, -0.7684, -0.7705, -0.7725, -0.7743, -0.7758, -0.7769, 
#             -0.7775, -0.7775]
# F0a_poly = nppoly.Polynomial.fit(p_nodes, F0a_data, 5)

# # For polynomial fit method Tc
# # From Regan, Wiman, Sauls arXiv:1908.04190 Table 2
# # Doesn't work at the moment, some unit miunserstanding or misprint?
# # a1 = [-9.849e-3, -5.043e-2, 2.205e-2, -2.557e-2, 5.023e-2 -2.769e-2]
# # a_list = [a1[::-1]] # polyfit wants highest power first
# # b1_poly = np.polynomial.Polynomial(a1)

# ###############################################################################################################
# ##########            Greywall, PLTS temperature scales polynomial coefficients                 ###############
# ###############################################################################################################


# # From Greywall 1986
# T_pcp_mK = 2.273
# p_pcp_bar = 21.22
# p_A_bar = 34.338
# # Greywall 1986 Tc polynomial fit cofficients and constructor
# aTc_G = [0.92938375, 0.13867188, -0.69302185e-2, 0.25685169e-3, -0.57248644e-5, 0.5301091e-7]
# Tc_poly_Greywall = nppoly.Polynomial(aTc_G)
# # Greywall 1986 TAB polynomial fit cofficients and constructor
# aTAB_G = [T_pcp_mK, -0.10322623e-1, -0.53633181e-2, 0.83437032e-3, -0.61709783e-4,  0.17038992e-5]
# TAB_poly_Greywall = nppoly.Polynomial(aTAB_G)

# # Melting curve poly coeffs (Greywall) 0.9 - 5.6 mK
# aPmelt = [-0.19632970e-1, 0.61880268e-1, -0.781103055e-1, 0.13050600, -0.43519381e-1, 
#           0.13752791e-6, -0.17180436e-6, -0.220939060e-9, -0.85450245e-12] 
# Pmelt_poly_Greywall = nppoly.Polynomial(aPmelt)

# # Greywall scale to PLTS scale, 6 order polynomial coefficients (Parpia et al 2022)
# GtoPLTS6 = [-0.14265343150487, 1.2810635032153, -0.22689947807354, 0.084337673002034, 
#             -0.016928990685839, 0.0017611612884063, -7.4461876859237e-5]
# GtoPLTS6_low_poly = np.polynomial.Polynomial(GtoPLTS6)

# # Greywall scale to PLTS scale, 9 order polynomial coefficients (Parpia et al 2022)
# GtoPLTS9 = [0.020353327019475, 0.96670033496024, 0.0019559314169033, -9.5551084662924e-5, 
#             3.2167457655106e-6, -7.0097586342143e-8, 9.6909878738352e-10,
#             -8.2126513949290e-12, 3.8886762300964e-14, -7.8713540127550e-17]
# GtoPLTS9_high_poly = np.polynomial.Polynomial(GtoPLTS9)


# # PLTS scale coefficients for polynomial method for Tc and T_AB, from Parpia et al 2022
# d_c = np.array([0.90972399274531, 0.14037182852625, -0.0074017331747577, 2.8617547367067e-4,-6.5064429600510e-6, 6.0754459040296e-8])
# c_AB = np.array([-26.864685876026, 5.2647866128370, -0.37617826876151, 0.013325635880953, -2.3510107585468e-4, 1.6519539175010e-6])

# Tc_poly_PLTS = nppoly.Polynomial(d_c)
# TAB_poly_PLTS = nppoly.Polynomial(c_AB)

# # Melting curve coefficients, PLTS, 0.9 - 5.6 mK, p in millibar
# b_melt = np.array([2.4492188707375, -0.026522783946809, -2.7665556176967e-5, -2.2800036357249e-7, 
#        3.9890194953355e-10, 1.2845930171276e-11, -6.9521369379387e-13, -1.5658128429388e-14, 
#        -1.1687750824147e-16, -3.0199721850282e-19] )
# Tmelt_poly_PLTS = nppoly.Polynomial(b_melt/1e3) # convert to bar

# ###############################################################################################################
# ##########            Order parameter useful basuis vectors, matrices and operations            ###############
# ###############################################################################################################
# # Basis vectors

# e = []
# e.append(np.array([0, 0, 1]))
# e.append(np.array([1, 1j, 0])/np.sqrt(2))
# e.append(np.array([1, -1j, 0])/np.sqrt(2))

# # Some basis matrices
# z3 = np.zeros((3,3))
# id3 = np.identity(3)
# D_A = np.outer(e[0], e[1])
# D_B = id3/np.sqrt(3)
# D_planar = (id3 - np.outer(e[0], e[0]))/np.sqrt(2)
# D_polar = np.outer(e[0], e[0])

# # Lowest barrier from A phase by exhaustive search
# D_low = np.array([[-0.16903589-0.2054976j,  -0.24395354-0.43379841j,  0.0228508 -0.06064158j],
#  [-0.03924275-0.003804j,    0.05325473-0.02309472j,  0.6362785 -0.39972627j],
#  [ 0.07959769-0.05774015j,  0.24372012-0.19001106j,  0.04900674-0.0131628j ]])

# # Antihermitian generators
# T_xy = np.array([[0,1,0],[-1,0,0],[0,0,0]])
# T_yz = np.array([[0,0,0],[0,0,1],[0,-1,0]])
# T_xz = np.array([[0,0,1],[0,0,0],[-1,0,0]])
# # pi/2 rotations
# O_xy = np.array([[0,1,0],[-1,0,0],[0,0,1]])
# O_yz = np.array([[1,0,0],[0,0,1],[0,-1,0]])
# O_xz = np.array([[0,0,1],[0,1,0],[-1,0,0]])

# # Dictionary of phases.
# inert_phases = ["B", "planar", "polar", "alpha", "bipolar", "A", "Ay", "Ayy", "Az", "beta", "gamma" ]

# R_arr_list = [np.array([1, 1, 1/3, 1/3, 1/3]),
#               np.array([1, 1, 1/2, 1/2, 1/2]),
#               np.array([1, 1, 1,   1,   1]),
#               np.array([0, 1, 1/3, 1/3, 1/3]),
#               np.array([0, 1, 1/2, 1/2, 1/2]),
#               np.array([0, 1, 0,   1,   1]),
#               np.array([0, 1, 0,   1,   1]),
#               np.array([0, 1, 0,   1,   1]),
#               np.array([0, 1, 0,   1,   1]),
#               np.array([0, 1, 1,   1,   0]),
#               np.array([0, 1, 0,   1,   0])]

# R_dict = dict(zip(inert_phases, R_arr_list))


# # Need to have a normalised OP class separate from phases

# D_dict = { "B"       : id3/np.sqrt(3),
#            "planar"  : (id3 - np.outer(e[0], e[0]))/np.sqrt(2),
#            "polar"   : np.matmul(np.matmul(O_xz, D_polar), np.transpose(O_xz) ), 
#            "alpha"   : np.diag(np.array([1, np.exp(1j*np.pi/3), np.exp(2j*np.pi/3)])/np.sqrt(3)),
#            "bipolar" : np.diag(np.array([1, 1j, 0])/np.sqrt(2)),
#            "A"       : np.matmul(O_xz, D_A),
#            "Ay"      : -1j*np.matmul(O_yz,np.matmul(D_A, O_xz)), #-1j*np.matmul(O_yz, D_A),
#            "Ayy"     : np.matmul(np.matmul(O_yz,D_A), O_yz), #-1j*np.matmul(O_yz, D_A),
#            "Az"      : D_A,
#            "beta"    : np.matmul(np.outer(e[1], e[0]), O_xz), 
#            "gamma"   : np.outer(e[1], e[1])
#            }




# Experimental data functions
def Tc_mK_expt(p):
    """
    Critical temperature in mK, using Greywall 1986 data, polynomial interpolated.
    T scale is set by current value of he3_tools.DEFAULT_T_SCALE.
    """
    # if scale == "PLTS":
    if DEFAULT_T_SCALE == "PLTS":

        return h3d.Tc_poly_PLTS(p)
    else:
        return h3d.Tc_poly_Greywall(p)

def Tc_mK(p):
    """ Wrapper for Tc_mK_expt(p).
    """
    # Tc_interp = np.interp(p, p_nodes, Tc_data_mK)
    # # if scale == "PLTS":
    # if DEFAULT_T_SCALE == "PLTS":
    #     return T_G_to_PLTS(Tc_interp)
    # else:
    #     return Tc_interp
    
    return Tc_mK_expt(p)


def T_mK(t, p):
    """Converts reduced temperature to temperature in mK.
    """
    return t * Tc_mK_expt(p)

def TAB_mK_expt(p):
    """
    AB equilibrium temperature in mK, using Greywall 1986 data, polynomial interpolated.
    T scale is set by current value of he3_tools.DEFAULT_T_SCALE.
    If pressure is less than polycritical point, he3_tools.p_pcp_bar, then np.nan
    is returned.
    """
    if DEFAULT_T_SCALE == "PLTS":
        TAB = h3d.TAB_poly_PLTS(p)
    else:
        TAB = h3d.TAB_poly_Greywall(p - h3d.p_pcp_bar)
    if isinstance(p, np.ndarray):
        TAB[p<h3d.p_pcp_bar] = np.nan
    else:
        if p < h3d.p_pcp_bar:
            TAB = np.nan
    return TAB

def tAB_expt(p):
    """
    AB equilibrium reduced temperature, using Greywall 1986 data, polynomial interpolated.
    If pressure is less than polycritical point, he3_tools.p_pcp_bar, then np.nan
    is returned.
    """
    return TAB_mK_expt(p)/Tc_mK_expt(p)

def p_melt(T_mK):
    # Melting pressure - needs updating.
    return h3d.p_A_bar * np.ones_like(T_mK)

def T_melt_PLTS(p):
    return h3d.Tmelt_poly_PLTS(p)

# Temperature scale convertor, Greywall to PLTS, ninth order polynomial 
def T_GtoPLTS(TG):  
    # return np.poly1d(np.flip(GtoPLTS9_high_poly.coef))(TG)
    return nppoly.Polynomial(h3d.GtoPLTS9_high_poly.coef)(TG)

def npart(p):
    """Particle density at pressure p.
    """
    # return np.interp(p, p_nodes, np_data)
    return np.interp(p, h3d.p_nodes, h3d.np_data)

def mstar_m(p):
    """Effective mass ratio at pressure p.
    """
    # return np.interp(p, p_nodes, mstar_m_data)
    return np.interp(p, h3d.p_nodes, h3d.mstar_m_data)

def vf(p):
    """Fermi velocity at pressure p.
    """
    # return np.interp(p, p_nodes, vf_data)
    return np.interp(p, h3d.p_nodes, h3d.vf_data)

def xi0(p):
    """Zero T Cooper pair correlation length at pressure p (nm).
    """
    # return np.interp(p, p_nodes, xi0_data)
    return np.interp(p, h3d.p_nodes, h3d.xi0_data)

def F0a(p):
    """Landau parameter $F_0^a$."""
    return h3d.F0a_poly(p)

def xi(t, p):
    """Ginzburg Landau correlation length at pressure p.
    """
    return h3c.xiGL_const*xi0(p)/(-alpha_norm(t))**0.5

def xi_delta(t, p):
    """BCS correlation length at pressure p.
    """
    return h3c.xiGL_const*xi0(p)/(-alpha_bcs(t))**0.5
    

def N0(p):
    """Density of states at Fermi surface, units nm^{-3} J^{-1}.
    """
    return npart(p) * 1.5 / (h3c.mhe3_kg * mstar_m(p) * vf(p)**2)

# Theory functions
def f_scale(p):

    """Free energy density units Joule per nm3 .
    """
    # return (1/3) * N0(p) * (2 * np.pi * kB * T_mK(1, p) * 1e-3)**2
    return (1/3) * N0(p) * (h3c.kB * T_mK(1, p) * 1e-3)**2
    
def delta_beta_norm(p, n):
    """
    Strong coupling corrections to GL free energy material parameters, in units of 
    the modulus of the first weak coupling parameter.
    
    Various calculations in the literature exist. The calculations and the 
    method for evaluating at a given pressure are set by the global variable 
    DEFAULT_SC_CORRS.

    'RWS19' [default]: Regan, Wiman and Sauls 2019 table I (said in that paper 
    to be taken from Wiman thesis) of strong coupling coefficients against presure 
    in bar (every 2 bar between 0 and 34). The default evaluation method is 
    interpolation.

    'RWS19-poly': data as above, fitted to 5th degree polynomials.
    
    'Choi-interp': Choi et. al 2007 table, interpolated. 

    'Choi-poly': Choi et. al 2007 table, fitted to polynomials of various degrees, 
    as specified in Wiman & Sauls 2015.
    
    'WS15': Choi et al data 2007, polynomial fit coefficients as given in WS15.
    
    'WS15-poly': same as WS15.

    'Wiman-thesis': From beta curves in Figure 8, grabbed using WPD, interpolated.
    """

    if DEFAULT_SC_CORRS in sc_corrs_poly:
        db = h3d.dbeta_data_dict[DEFAULT_SC_CORRS][n-1](p)
    else:
        db = delta_beta_norm_interp(p, n)
    
    if DEFAULT_SC_ADJUST:
            db *= np.exp(-sc_adjust_fun(p))
    return db

def delta_beta_norm_asarray(p):
    """
    Strong coupling corrections to material parameters, in units of 
    the modulus of the first weak coupling parameter, supplied as a (1,5) array.
    """
    delta_beta_norm_list = [ delta_beta_norm(p, n) for n in range(1,6)]
    return np.array(delta_beta_norm_list)

def delta_beta_norm_interp(p, n): 
    """Interpolation methods for strong coupling corrections.
    """
    if DEFAULT_SC_CORRS in sc_corrs_interp:
        p_nodes_beta = h3d.dbeta_data_dict[DEFAULT_SC_CORRS][0]
        c_list = h3d.dbeta_data_dict[DEFAULT_SC_CORRS][1]
        return np.interp(p, p_nodes_beta, c_list[n-1])
    else:
        raise ValueError("No interpolation data for this value of DEFAULT_SC_CORRS:",
                         DEFAULT_SC_CORRS)
        return

def delta_beta_norm_polyfit(p, n): 
    """Polynomial methods for strong couping corrections. 
    """
    if DEFAULT_SC_CORRS in sc_corrs_poly:
        return h3d.dbeta_data_dict[DEFAULT_SC_CORRS][n-1](p)
    else:
        raise ValueError("No polynomial for this value of DEFAULT_SC_CORRS", 
                         DEFAULT_SC_CORRS)
        return

def alpha_bcs(t):
    """Fit to function giving BCS gap, asymptotes to (t - 1) asa t -> 1.
    """
    return (t**h3c.a_bcs - 1 )/h3c.a_bcs

def alpha_norm(t):
    """Quadratic material parameter
    """
    if DEFAULT_ALPHA_TYPE == "GL":
        return t - 1
    elif DEFAULT_ALPHA_TYPE == "BCS":
        return alpha_bcs(t)
    else:
        raise ValueError("DEFAULT_ALPHA_TYPE must be GL or BCS")
        return

def beta_norm(t, p, n):
    """Complete material parameter including strong coupling correction, within units of 
    f_scale/(2 * np.pi * kB * Tc)**2
    """ 
    b = beta_norm_wc(n)
    # else:
    #     raise ValueError("beta_norm: n must be between 1 and 5")
    return h3c.beta_const*(b + t * delta_beta_norm(p, n))


def beta_norm_asarray(t, p):
    beta_norm_list = [ beta_norm(t, p, n) for n in range(1,6)]
    return np.array(beta_norm_list)

def mat_pars(t,p):
    """All 6 bulk material parameters alpha, beta_a as a numpy array.
    """
    pars = np.zeros((6,))
    pars[0] = alpha_norm(t)
    pars[1:] = beta_norm_asarray(t, p)
    return pars

def beta_phase_norm(t, p, phase):
    """Effective single beta parameter in a given phase.
    """
    return np.sum( beta_norm_asarray(t, p) * h3b.R_dict[phase] )

def f_phase_norm(t, p, phase):
    """Normalised free energy in a givenn phase.
    """
    return -0.25* alpha_norm(t)**2 /( beta_phase_norm(t, p, phase))

def delta_phase_norm(t, p, phase):
    """Normalised gap parameter in a given phase.
    """
    return np.sqrt(- alpha_norm(t)/(2 * beta_phase_norm(t, p, phase)))

def delta_wc(t):
    return np.sqrt(- alpha_norm(t)/(2 * h3c.beta_const))

def beta_A_norm(t, p):    
    """Normalised effective single beta parameter for A phase.
    """
    return beta_norm(t, p, 2) +  beta_norm(t, p, 4) + beta_norm(t, p, 5)

def beta_B_norm(t, p):
    """Normalised effective single beta parameter for B phase.
    """
    return beta_norm(t, p, 1) + beta_norm(t, p, 2) + (beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/3


def beta_planar_norm(t, p):
    """Normalised effective single beta parameter for planar phase.
    """
    return beta_norm(t, p, 1) + beta_norm(t, p, 2) + (beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/2

def beta_polar_norm(t, p):
    """Normalised effective single beta parameter for polar phase.
    """
    return beta_norm(t, p, 1) + beta_norm(t, p, 2) + (beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))

def f_A_norm(t, p):
    """Normalised free energy density for A phase, in units of (1/3) N(0) (k_B T_c)^2.
    """
    return -0.25* alpha_norm(t)**2 /( beta_A_norm(t, p))
    
def f_B_norm(t, p):
    """Normalised free energy density for B phase, in units of (1/3) N(0) (k_B T_c)^2.
    """
    return -0.25* alpha_norm(t)**2 /( beta_B_norm(t, p))

def f_planar_norm(t, p):
    """Normalised free energy density for planar phase, in units of (1/3) N(0) (k_B T_c)^2.
    """
    return -0.25* alpha_norm(t)**2 /( beta_planar_norm(t, p))

def delta_A_norm(t, p):
    """Gap parameter for A phase, normalised to (kB * Tc)
    """    
    return np.sqrt(- alpha_norm(t)/(2 * beta_A_norm(t,p)))


def delta_B_norm(t, p):
    """Gap parameter for B phase, normalised to (kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_B_norm(t, p)))

def delta_planar_norm(t, p):
    """Gap parameter for planar phase, normalised to (kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_planar_norm(t, p)))

def delta_polar_norm(t, p):
    """Gap parameter for planar phase, normalised to (kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_polar_norm(t, p)))

def t_AB(p):
    """ AB transition temperature at pressure p, normalised to Tc.
    """
    t_ab_val = (1/3)/ (delta_beta_norm(p, 1) 
                       + (delta_beta_norm(p, 3) 
                       - 2*delta_beta_norm(p, 4) 
                       - 2*delta_beta_norm(p, 5))/3) 
    
    if isinstance(t_ab_val, np.ndarray):
        t_ab_val[t_ab_val > 1] = np.nan
    elif t_ab_val > 1:
        t_ab_val = np.nan
            
    return  t_ab_val

def tAB(p):
    """Synonym for t_AB.
    """
    return t_AB(p)

def TAB_mK(p):
    """ AB transition temperature at pressure p, in mK
    """
    return tAB(p) * Tc_mK_expt(p)

# Generate SC adjustment factor
def logf_poly():
    p = np.linspace(h3d.p_pcp_bar, 34, 100)
    global DEFAULT_SC_ADJUST
    tmp = DEFAULT_SC_ADJUST
    DEFAULT_SC_ADJUST = False
    logf = np.log(tAB_expt(p)/t_AB(p))
    DEFAULT_SC_ADJUST=tmp
    return nppoly.Polynomial.fit(p, logf, 2)

def sc_adjust_fun(p):
    sc_corr_adj_pol = logf_poly()
    adj_exp = sc_corr_adj_pol(p)
    if DEFAULT_SC_ADJUST_MIN_P=="p_pcp_bar":
        p0 = h3d.p_pcp_bar
    else:
        p0 = 0
    if isinstance(adj_exp, np.ndarray) or isinstance(adj_exp, pd.Series):
        adj_exp[p < p0] = 0.0
    else: # assume float
        if p < p0:
            adj_exp = 0.0
    return adj_exp

def mass_B_norm(t, p, JC):
    """B phase masses for mode with spin parity JC
    """
    bb = beta_B_norm(t, p)
    
    if JC == "1-":        
        m2 = (- beta_norm(t, p, 1)+ (beta_norm(t, p, 4) - beta_norm(t, p, 3) - beta_norm(t, p, 5))/3) / bb
    elif JC == "2+":
        m2 = ((beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/3) / bb
    elif JC == "2-":
        m2 = (- beta_norm(t, p, 1)) / bb

    return np.sqrt(m2)        


def critical_radius(t, p, sigma=0.95, dim=3):
    """Radius of critical bubble, in nm.  
    Ideally will optionally use function to get 
    surface tension. Uses approximation."""
    # if isinstance(sigma_fun, float):
    sigma_AB = sigma*np.abs(f_B_norm(t,p))*xi(t,p)
    # elif isinstance(sigma_fun, np.ndarray):
        # sigma_AB = sigma_fun*np.abs(f_B_norm(t,p))*xi(t,p)
    
    return (dim-1)*sigma_AB/np.abs(f_A_norm(t,p) - f_B_norm(t,p))


# def tr(A):
#     """
#     Trace of order parameter array on last two indices

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     d = A.ndim    
#     return np.trace(A, axis1=d-2,  axis2=d-1 )


# def trans(A):
#     """
#     Transpose of order parameter array on last two indices

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     d = A.ndim    
#     return np.swapaxes(A, d-2,  d-1 )


# def hconj(A):
#     """
#     Hermitian conjugate of order parameter array on last two indices

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     return trans(A).conj()


# def mult_hconj(A):
#     """
#     Multiplies A by its Hermitian conjugate on last two indices

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     ndarray
#         Same shape as input array.

#     """
#     return np.matmul(A, hconj(A))

# def mult_trans(A):
#     """
#     Multiplies A by its transpose on last two indices

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     ndarray
#         Same shape as input array.

#     """
#     return np.matmul(A, trans(A))

# def mat2vec(A):
#     sz = A.size 
#     return A.reshape(sz)

# def vec2mat(Av):
#     sz = Av.size 
#     return Av.reshape((sz//9, 3, 3))


# def eig_vals(A):
#     """
#     Extracts eigenvalues of array A

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     ndarray
#         shape (m,n,...,3).

#     """
#     A_dim = A.ndim
#     A_shape = A.shape
#     A_size = A.size
#     A_flatter = A.reshape(A_size//9, 3, 3)
#     e_vals = np.zeros((A_size//9, 3), dtype=complex)
#     for n, A_el in enumerate(A_flatter):
#         e_vals[n,:] = sl.eigvals(A_el)
#     return e_vals.reshape(A_shape[0:A_dim-2] + (3,))


# def eig_orbital(A):
#     """
#     Extracts eigenvalues of array A^dagger A, relevant for angular momentum vectors

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     ndarray
#         shape (m,n,...,3).

#     """
#     H = np.matmul(hconj(A), A)
#     H_dim = H.ndim
#     H_shape = H.shape
#     H_size = H.size
#     H_flatter = H.reshape(H_size//9, 3, 3)
#     e_vals = np.zeros((H_size//9, 3), dtype=float)
#     for n, H_el in enumerate(H_flatter):
#         e_vals[n,:] = sl.eigvals(H_el)
#     e_vals.sort()
#     return e_vals.reshape(H_shape[0:H_dim-2] + (3,))


# def eig_spin(A):
#     """
#     Extracts eigenvalues of array A A^dagger, relevant for spin vectors

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter array.

#     Returns
#     -------
#     ndarray
#         shape (m,n,...,3).

#     """
#     H = np.matmul(A, hconj(A))
#     H_dim = H.ndim
#     H_shape = H.shape
#     H_size = H.size
#     H_flatter = H.reshape(H_size//9, 3, 3)
#     e_vals = np.zeros((H_size//9, 3), dtype=float)
#     for n, H_el in enumerate(H_flatter):
#         e_vals[n,:] = sl.eigvals(H_el)
#     e_vals.sort()
#     return e_vals.reshape(H_shape[0:H_dim-2] + (3,))


# def args_parse(*args):
#     n_el = len(args)
#     if n_el == 2:
#         if isinstance(args[1], np.ndarray):
#             t = None
#             p = None
#             al = args[0]
#             bn = args[1]
#         else:
#             t, p = args
#             al = alpha_norm(t)
#             bn = beta_norm_asarray(t, p)
#     else:
#         t = None
#         p = None
#         al = args[0][0]
#         bn = args[0][1:]
#     return al, bn, t, p


# # def U(A, alpha_norm, beta_norm_arr):
# def U(A, *args):
#     """
#     Bulk free energy for superfluid He3

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter.
#     alpha_norm : float
#         alpha parameter, normalised.
#     beta_norm_arr : ndarray, shape (5,1)
#         beta parameters, normalised.

#     Returns
#     -------
#     float or ndarray
#         Normalised bulk free energy.

#     """
#     al, bn, _, _ = args_parse(*args)

#     dim = A.ndim
    
#     A_T = np.swapaxes(A, dim-2, dim-1)
#     A_C = np.conj(A)
#     A_H = np.conj(A_T)

#     Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
#     Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
#     Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
#     Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
#     Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
#     Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
#     return (Un0 + Un1 + Un2 + Un3 + Un4 + Un5).real


# # def U(A, alpha_norm, beta_norm_arr):
# def U_terms(A, *args):
#     """
#     Bulk free energy for superfluid He3

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...,3,3)
#         order parameter.
#     alpha_norm : float
#         alpha parameter, normalised.
#     beta_norm_arr : ndarray, shape (5,1)
#         beta parameters, normalised.

#     Returns
#     -------
#     float or ndarray
#         Normalised bulk free energy.

#     """
#     al, bn, _, _ = args_parse(*args)

#     dim = A.ndim
    
#     A_T = np.swapaxes(A, dim-2, dim-1)
#     A_C = np.conj(A)
#     A_H = np.conj(A_T)

#     Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
#     Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
#     Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
#     Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
#     Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
#     Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
#     return np.array([Un0, Un1, Un2, Un3, Un4, Un5])


# def dU_dA(A, *args):
#     """
#     Derivative of bulk free energy for superfluid He3.')

#     Parameters
#     ----------
#     A : ndarray dtype complex, shape (m,n,...3,3)
#         order parameter.
#     alpha_norm : float
#         alpha material parameter, dimensionless.
#     beta_norm_arr : ndarray, shape (5,1)
#         beta material parameters, dimensionless.

#     Returns
#     -------
#     float or ndarray
#         Normalised derivative of bulk free energy.

#     """
#     al, bn, _, _ = args_parse(*args)

#     dim = A.ndim
    
#     A_T = np.swapaxes(A, dim-2, dim-1)
#     A_C = np.conj(A)
#     A_H = np.conj(A_T)

#     dUn0 = al * A
        
#     dUn1 = 2 * bn[0] *  np.matmul(np.multiply.outer(h3m.tr( np.matmul(A , A_T)), h3b.id3 ) , A_C)
#     dUn2 = 2 * bn[1] *  np.matmul(np.multiply.outer(h3m.tr( np.matmul(A , A_H)), h3b.id3 ) , A )
#     dUn3 = 2 * bn[2] *  np.matmul(A , np.matmul(A_T , A_C))
#     dUn4 = 2 * bn[3] *  np.matmul(A , np.matmul(A_H   , A))
#     dUn5 = 2 * bn[4] *  np.matmul(A_C , np.matmul(A_T , A))
#     return dUn0 + dUn1 + dUn2 + dUn3 + dUn4 + dUn5


# def line_section(X, D, t, p, scale=None, n=500):
#     """
#     Retuens a line section in order parameter space, along with the free energy 
#     along it, and the parameter of the line.

#     Parameters
#     ----------
#     X : string or Complex array shape (3,3)
#         Start point of line. If a string, specifies the inert phase.
#     D : string or Complex array shape (3,3)
#         If a string, the end point inert phase. If array, the direction of the line
#     t : float
#         Reduced temperature in terms of Tc(p).
#     p : float
#         Pressure in bar.
#     scale : float, optional
#         Multiplies direction matrix. The default is None, in which case the 
#         scale is the weak coupling gap scale delta_wc(t).
#     n : integer, optional
#         Number of points on the line. The default is 500.

#     Returns
#     -------
#     v : numpy.ndarray shape (n,)
#         Line parameter, between 0 and v_max=1.5.
#     A_XD : Complex array shape (n,3,3)
#         Order parameter along line.
#     U_XD : float array shape (n,)
#         Free energy along line.

#     """
#     if scale is None:
#         scale = delta_wc(t)
#     v_max = 1.5
#     v = np.linspace(0,1,n)*v_max
  
#     if isinstance(X, str):
#         delta_X = delta_phase_norm(t, p, X)
#         X = h3b.D_dict[X]
#         if isinstance(D, str):
#             D = h3b.D_dict[D] - X
#             D = D/h3m.norm(D)

#     A_XD = np.multiply.outer( np.ones_like(v) , X)*delta_X + np.multiply.outer( v , D)*scale
#     U_XD = U(A_XD, (t-1), beta_norm_asarray(t, p) )
    
#     return v, A_XD, U_XD

# def norm(D, n=1):
#     """
#     n-Norm of complex square array D, defined as $tr(D D^\dagger)^{n/2}$

#     Parameters
#     ----------
#     D : Complex array shape 
#         Order parameter array.

#     Returns
#     -------
#     float
#         n-Norm of array.

#     """
#     # D_H = np.conj(np.swapaxes(D, dim-2, dim-1))
#     D_H = h3m.hconj(D)
#     norm2 = h3m.tr(np.matmul(D, D_H)).real
#     if n == 1:
#         return norm2
#     else:
#         return norm2**(n/2)


# def inner(X, Y):
#     """
#     Inner product of complex square arrays X, Y

#     Parameters
#     ----------
#     X,Y : Complex array shape 
#         Order parameter array.

#     Returns
#     -------
#     float, dtype complex
#         Inner product.

#     """
#     dim = Y.ndim
#     # X_H = np.conj(np.swapaxes(Y, dim-2, dim-1))
#     Y_H = np.conj(np.swapaxes(Y, dim-2, dim-1))
#     return (h3m.tr(np.matmul(X, Y_H))).real


# def distance(X, Y):
#     """
#     Distance between (3,3) complex arrays defined by norm.    

#     Parameters
#     ----------
#     X : Complex array shape (3,3)
#         Order parameter array.
#     Y : Complex array shape (3,3)
#         Order parameter array.

#     Returns
#     -------
#     float
#         Distance between X and Y.

#     """
#     D = X - Y
#     return norm(D)


# def project(X, Y):
#     """
#     Projects X onto Y, X, Y (3,3) complex arrays.
#     X - (X,Y) Y/(Y,Y)
       

#     Parameters
#     ----------
#     X : Complex array shape (3,3)
#         Order parameter array.
#     Y : Complex array shape (3,3)
#         Order parameter array.

#     Returns
#     -------
#     float
#         X projected onto Y.

#     """

#     return X -  np.multiply.outer(inner(X,Y)/inner(Y,Y), Y )


# def lengths(X, Y):
#     """
#     Length of X along Y and orthoginal to Y
#     Projects X onto Y, X, Y (3,3) complex arrays.
#     X - (X,Y) Y/(Y,Y)
       

#     Parameters
#     ----------
#     X : Complex array shape (n,m,..3,3)
#         Order parameter array.
#     Y : Complex array shape (3,3)
#         Order parameter array.

#     Returns
#     -------
#     float
#         X projected onto Y.

#     """
#     print((inner(X,Y)/inner(Y,Y)).shape)
#     print(Y.shape)
#     X_Y = X - np.multiply.outer(inner(X,Y)/inner(Y,Y), Y ) 

#     return np.array([norm(X_Y), norm(X - X_Y)]).T


# for x in ["DEFAULT_SC_ADJUST", "DEFAULT_SC_CORRS", "DEFAULT_T_SCALE"]:
#     report_setting(x)


    
