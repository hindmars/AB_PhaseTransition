#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:04:39 2021

@author: hindmars
"""

import numpy as np
import scipy.special as sp

import scipy.constants as c


cphy = c.physical_constants

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


Tc_data_mK = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 
      2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.48]

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


# For polynomial fit method
# From Regan, Wiman, Sauls arXiv:1908.04190 Table 2
# Doesn't work at the moment, some unit miunserstanding or misprint?
a1 = [9.849e-3, -5.043e-2, 2.205e-2, -2.557e-2, 5.023e-2 -2.769e-2]
a_list = [a1[::-1]] # polyfit wants highest power first


zeta3 = sp.zeta(3)
# beta_const = 7 * zeta3/(240 * np.pi**2)
beta_const = 7 * zeta3/(80 * np.pi**2)
xiGL_const = np.sqrt(7 * zeta3 /    20)

# Helium 3 mass in u
mhe3_u = 3.0160293 
mhe3_kg = mhe3_u * cphy["atomic mass constant"][0]

kB = c.k
hbar = c.hbar

# Order parameter directions

e = []
e.append(np.array([0, 0, 1]))
e.append(np.array([1, 1j, 0])/np.sqrt(2))
e.append(np.array([1, -1j, 0])/np.sqrt(2))

id3 = np.identity(3)
D_A = np.outer(e[0], e[1])
D_B = id3/np.sqrt(3)
D_planar = (id3 - np.outer(e[0], e[0]))/np.sqrt(2)
D_polar = np.outer(e[0], e[0])

# Lowest barrier by exhaustive search
D_low = np.array([[-0.16903589-0.2054976j,  -0.24395354-0.43379841j,  0.0228508 -0.06064158j],
 [-0.03924275-0.003804j,    0.05325473-0.02309472j,  0.6362785 -0.39972627j],
 [ 0.07959769-0.05774015j,  0.24372012-0.19001106j,  0.04900674-0.0131628j ]])

# Dictionary of phases

O_xz = np.array([[0,0,1],[0,1,0],[-1,0,0]])

inert_phases = ["B", "planar", "polar", "alpha", "bipolar", "A", "beta", "gamma" ]

R_arr_list = [np.array([1, 1, 1/3, 1/3, 1/3]),
              np.array([1, 1, 1/2, 1/2, 1/2]),
              np.array([1, 1, 1,   1,   1]),
              np.array([0, 1, 1/3, 1/3, 1/3]),
              np.array([0, 1, 1/2, 1/2, 1/2]),
              np.array([0, 1, 0,   1,   1]),
              np.array([0, 1, 1,   1,   0]),
              np.array([0, 1, 0,   1,   0])]

R_dict = dict(zip(inert_phases, R_arr_list))

D_dict = { "B"       : id3/np.sqrt(3),
           "planar"  : (id3 - np.outer(e[0], e[0]))/np.sqrt(2),
           "polar"   : np.matmul(np.matmul(O_xz, D_polar), np.transpose(O_xz) ), 
           "alpha"   : np.diag(np.array([1, np.exp(1j*np.pi/3), np.exp(2j*np.pi/3)])/np.sqrt(3)),
           "bipolar" : np.diag(np.array([1, 1j, 0])/np.sqrt(2)),
           "A"       : np.matmul(O_xz, D_A),
           "beta"    : np.matmul(np.outer(e[1], e[0]), O_xz), 
           "gamma"   : np.outer(e[1], e[1])
           }


# Functions
def Tc_mK(p):
    return np.interp(p, p_nodes, Tc_data_mK)


def T_mK(t, p):
    """Converts reduced temperature to temperature in mK.
    """
    return t * Tc_mK(p)


def npart(p):
    """Particle density at pressure p.
    """
    return np.interp(p, p_nodes, np_data)


def mstar_m(p):
    """Effective mass ratio at pressure p.
    """
    return np.interp(p, p_nodes, mstar_m_data)


def vf(p):
    """Fermi velocity at pressure p.
    """
    return np.interp(p, p_nodes, vf_data)


def xi0(p):
    """Zero T Cooper pair correlation length at pressure p.
    """
    return np.interp(p, p_nodes, xi0_data)

def xi(t, p):
    """Ginzburg Landau correlation length at pressure p.
    """
    return xiGL_const*xi0(p)/(1-t)**0.5

def N0(p):
    """Density of states at Fermi surface, units nm^{-3} J^{-1}.
    """
    return npart(p) * 1.5 / (mhe3_kg * mstar_m(p) * vf(p)**2)

def f_scale(p):
    """Free energy density units Joule per nm3 (?).
    """
    # return (1/3) * N0(p) * (2 * np.pi * kB * T_mK(1, p) * 1e-3)**2
    return (1/3) * N0(p) * (kB * T_mK(1, p) * 1e-3)**2
    
def delta_beta_norm(p, n, method="interp"):
    """Strong coupling corrections to material parameters, in units of 
    the modulus of the first weak coupling parameter.
    """
    if method == "interp":
        return delta_beta_norm_interp(p, n)
    elif method == "polyfit":
        return delta_beta_norm_polyfit(p, n)
    else:
        raise ValueError("error: strong coupling parameter method must be interp or polyfit")
        
    return

def delta_beta_norm_interp(p, n): 
    """Interpolation method.
    """
    return np.interp(p, p_nodes, c_list[n-1])


def delta_beta_norm_polyfit(p, n): 
    """Polynomial method, not implemented yet.
    """
    if n==1:
        return np.poly1d(a_list[n-1], len(a_list[n-1]))(p)
    else:
        raise ValueError("only beta_1 implemented as yet")
        return

def alpha_norm(t):
    return t - 1

def beta_norm(t, p, n):
    """Complete material parameter including strong coupling correction, in units of 
    f_scale/(2 * np.pi * kB * T)**2
    """ 
    if n==1:
        b = -1
    elif 1 < n < 5:
        b = 2
    elif n==5:
        b = -2
    else:
        raise ValueError("beta_norm: n must be between 1 and 5")
    return beta_const*(b + t * delta_beta_norm(p, n))

def beta_norm_asarray(t, p):
    beta_norm_list = [ beta_norm(t, p, n) for n in range(1,6)]
    return np.array(beta_norm_list)

def beta_phase_norm(t, p, phase):
    return np.sum( beta_norm_asarray(t, p) * R_dict[phase] )

def f_phase_norm(t, p, phase):
    return -0.25* alpha_norm(t)**2 /( beta_phase_norm(t, p, phase))

def delta_phase_norm(t, p, phase):
    return np.sqrt(- alpha_norm(t)/(2 * beta_phase_norm(t, p, phase)))

def delta_wc(t):
    return np.sqrt(- alpha_norm(t)/(2 * beta_const))

def beta_A_norm(t, p):    
    """Material parameter for A phase.
    """
    return beta_norm(t, p, 2) +  beta_norm(t, p, 4) + beta_norm(t, p, 5)

def beta_B_norm(t, p):
    """Material parameter for B phase.
    """
    return beta_norm(t, p, 1) + beta_norm(t, p, 2) + (beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/3


def beta_planar_norm(t, p):
    """Material parameter for planar phase.
    """
    return beta_norm(t, p, 1) + beta_norm(t, p, 2) + (beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/2

def beta_polar_norm(t, p):
    """Material parameter for polar phase.
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
    """Gap parameter for A phase, normalised to (2 * pi * kB * Tc)
    """    
    return np.sqrt(- alpha_norm(t)/(2 * beta_A_norm(t,p)))


def delta_B_norm(t, p):
    """Gap parameter for B phase, normalised to (2 * pi * kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_B_norm(t, p)))

def delta_planar_norm(t, p):
    """Gap parameter for planar phase, normalised to (2 * pi * kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_planar_norm(t, p)))

def delta_polar_norm(t, p):
    """Gap parameter for planar phase, normalised to (2 * pi * kB * Tc)
    """
    return  np.sqrt(- alpha_norm(t)/(2 * beta_polar_norm(t, p)))

def t_AB(p):
    """ AB transition temoerature at pressure p, normalised to Tc.
    """
    t_ab_val = (1/3)/ (delta_beta_norm(p,1) + (delta_beta_norm(p,3) - 2*delta_beta_norm(p,4) - 2*delta_beta_norm(p,5))/3) 
    
    if isinstance(t_ab_val, np.ndarray):
        t_ab_val[t_ab_val > 1] = np.nan
    elif t_ab_val > 1:
        t_ab_val = np.nan
            
    
    return  t_ab_val


def mass_B_norm(t, p, JC):
    
    bb = beta_B_norm(t, p)
    
    if JC == "1-":        
        m2 = (- beta_norm(t, p, 1)+ (beta_norm(t, p, 4) - beta_norm(t, p, 3) - beta_norm(t, p, 5))/3) / bb
    elif JC == "2+":
        m2 = ((beta_norm(t, p, 3) + beta_norm(t, p, 4) + beta_norm(t, p, 5))/3) / bb
    elif JC == "2-":
        m2 = (- beta_norm(t, p, 1)) / bb

    return np.sqrt(m2)        


def tr(a):
    d = a.ndim    
    return np.trace(a, axis1=d-2,  axis2=d-1 )


def U(A, alpha_norm, beta_norm_arr):
    """
    Bulk free energy for superfluid He3

    Parameters
    ----------
    A : ndarray
        order parameter.
    alpha_norm : float
        alpha parameter, normalised.
    beta_norm_arr : ndarray, shape (5,1)
        beta parameters, normalised.

    Returns
    -------
    float or ndarray
        Normalised bulk free energy.

    """
    bn = beta_norm_arr
    dim = A.ndim
    
    A_T = np.swapaxes(A, dim-2, dim-1)
    A_C = np.conj(A)
    A_H = np.conj(A_T)

    Un0 = alpha_norm * tr( np.matmul(A , A_H) )
    
    Un1 = bn[0] *  tr( np.matmul(A , A_T)) * tr(np.matmul(A_C , A_H) ) 
    Un2 = bn[1] *  tr( np.matmul(A , A_H) )**2
    Un3 = bn[2] *  tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
    Un4 = bn[3] *  tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
    Un5 = bn[4] *  tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
    return (Un0 + Un1 + Un2 + Un3 + Un4 + Un5).real


def dU_dA(A, alpha_norm, beta_norm_arr):
    """
    Derivative of bulk free energy for superfluid He3.')

    Parameters
    ----------
    A : ndarray
        order parameter.
    alpha_norm : float
        alpha parameter, normalised.
    beta_norm_arr : ndarray, shape (5,1)
        beta parameters, normalised.

    Returns
    -------
    float or ndarray
        Normalised derivative of bulk free energy.

    """
    bn = beta_norm_arr
    dim = A.ndim
    
    A_T = np.swapaxes(A, dim-2, dim-1)
    A_C = np.conj(A)
    A_H = np.conj(A_T)

    dUn0 = alpha_norm * A
        
    dUn1 = 2 * bn[0] *  np.matmul(np.multiply.outer(tr( np.matmul(A , A_T)), id3 ) , A_C)
    dUn2 = 2 * bn[1] *  np.matmul(np.multiply.outer(tr( np.matmul(A , A_H)), id3 ) , A )
    dUn3 = 2 * bn[2] *  np.matmul(A , np.matmul(A_T , A_C))
    dUn4 = 2 * bn[3] *  np.matmul(A , np.matmul(A_H   , A))
    dUn5 = 2 * bn[4] *  np.matmul(A_C , np.matmul(A_T , A))
    return dUn0 + dUn1 + dUn2 + dUn3 + dUn4 + dUn5


def line_section(X, D, t, p, scale=None, n=500):
    
    if scale is None:
        scale = delta_wc(t)
        
    v = np.linspace(0,1,n)*1.5
  
    if isinstance(X, str):
        delta_X = delta_phase_norm(t, p, X)
        X = D_dict[X]
        if isinstance(D, str):
            D = D_dict[D] - X

    A_XD = np.multiply.outer( np.ones_like(v) , X)*delta_X + np.multiply.outer( v , D)*scale
    U_XD = U(A_XD, (t-1), beta_norm_asarray(t, p) )
    
    return v, A_XD, U_XD

def norm(D):
    dim = D.ndim
    D_H = np.conj(np.swapaxes(D, dim-2, dim-1))
    return np.sqrt(tr(np.matmul(D, D_H))).real
    
def distance(X, Y):
    D = X - Y
    return norm(D)
    