#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 17:00:45 2022

Basis vectors and matrices for superfluid He3 order parameter.

@author: hindmars
"""

import numpy as np

###############################################################################################################
##########            Order parameter useful basuis vectors, matrices and operations            ###############
###############################################################################################################
# Basis vectors

e = []
e.append(np.array([0, 0, 1]))
e.append(np.array([1, 1j, 0])/np.sqrt(2))
e.append(np.array([1, -1j, 0])/np.sqrt(2))

# Some basis matrices
z3 = np.zeros((3,3))
id3 = np.identity(3)
D_A = np.outer(e[0], e[1])
D_B = id3/np.sqrt(3)
D_planar = (id3 - np.outer(e[0], e[0]))/np.sqrt(2)
D_polar = np.outer(e[0], e[0])


# Lowest barrier from A phase by exhaustive search
D_low = np.array([[-0.16903589-0.2054976j,  -0.24395354-0.43379841j,  0.0228508 -0.06064158j],
 [-0.03924275-0.003804j,    0.05325473-0.02309472j,  0.6362785 -0.39972627j],
 [ 0.07959769-0.05774015j,  0.24372012-0.19001106j,  0.04900674-0.0131628j ]])

# Antihermitian generators
T_xy = np.array([[0,1,0],[-1,0,0],[0,0,0]])
T_yz = np.array([[0,0,0],[0,0,1],[0,-1,0]])
T_xz = np.array([[0,0,1],[0,0,0],[-1,0,0]])
# pi/2 rotations
O_xy = np.array([[0,1,0],[-1,0,0],[0,0,1]])
O_yz = np.array([[1,0,0],[0,0,1],[0,-1,0]])
O_xz = np.array([[0,0,1],[0,1,0],[-1,0,0]])

# Dictionary of phases.
inert_phases = ["B", "planar", "polar", "alpha", "bipolar", "A", "Ay", "Ayy", "Az", "beta", "gamma" ]

R_arr_list = [np.array([1, 1, 1/3, 1/3, 1/3]),
              np.array([1, 1, 1/2, 1/2, 1/2]),
              np.array([1, 1, 1,   1,   1]),
              np.array([0, 1, 1/3, 1/3, 1/3]),
              np.array([0, 1, 1/2, 1/2, 1/2]),
              np.array([0, 1, 0,   1,   1]),
              np.array([0, 1, 0,   1,   1]),
              np.array([0, 1, 0,   1,   1]),
              np.array([0, 1, 0,   1,   1]),
              np.array([0, 1, 1,   1,   0]),
              np.array([0, 1, 0,   1,   0])]

R_dict = dict(zip(inert_phases, R_arr_list))


# Need to have a normalised OP class separate from phases

D_dict = { "B"       : id3/np.sqrt(3),
           "planar"  : (id3 - np.outer(e[0], e[0]))/np.sqrt(2),
           "polar"   : np.matmul(np.matmul(O_xz, D_polar), np.transpose(O_xz) ), 
           "alpha"   : np.diag(np.array([1, np.exp(1j*np.pi/3), np.exp(2j*np.pi/3)])/np.sqrt(3)),
           "bipolar" : np.diag(np.array([1, 1j, 0])/np.sqrt(2)),
           "A"       : np.matmul(O_xz, D_A),
           "Ay"      : -1j*np.matmul(O_yz,np.matmul(D_A, O_xz)), #-1j*np.matmul(O_yz, D_A),
           "Ayy"     : np.matmul(np.matmul(O_yz,D_A), O_yz), #-1j*np.matmul(O_yz, D_A),
           "Az"      : D_A,
           "beta"    : np.matmul(np.outer(e[1], e[0]), O_xz), 
           "gamma"   : np.outer(e[1], e[1]), 
           "B1"      : np.diag([-1, 1, 1])/np.sqrt(3), 
           "B2"      : np.diag([1, -1, 1])/np.sqrt(3), 
           "B3"      : np.diag([1, 1, -1])/np.sqrt(3), 
           "B4"      : np.diag([-1, -1, 1])/np.sqrt(3), 
           "B5"      : np.diag([-1, 1, -1])/np.sqrt(3), 
           "B6"      : np.diag([1, -1, -1])/np.sqrt(3), 
           "B7"      : np.diag([-1, -1, -1])/np.sqrt(3)
           }


