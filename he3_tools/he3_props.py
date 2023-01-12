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
    """
    Melting pressure , Greywall 1986 (equation A1).
    Needs PLTS to Greywall temperature converter
    """
    T = T_mK.copy()
    if DEFAULT_T_SCALE == "PLTS":
        T = T_PLTStoG(T_mK)
    # return h3d.p_A_bar * np.ones_like(T_mK)
    return h3d.Pmelt_poly_Greywall(T) + h3d.p_A_bar

# def T_melt_PLTS(p):
#     """This makes no sense - melting pressure is close to 34 bar. What is range of p?"""
#     return h3d.Tmelt_poly_PLTS(p)

# Temperature scale convertor, Greywall to PLTS, ninth order polynomial 
def T_GtoPLTS(TG):  
    return h3d.GtoPLTS9_high_poly(TG)

# Temperature scale convertor, invert Greywall to PLTS, ninth order polynomial 
def T_PLTStoG(TPLTS):  
    return h3d.PLTStoG9_high_poly(TPLTS)

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

def gH(p):
    """
    Material parameter for Zeeman energy, quadratic in magnetic field H, in 
    dimensionless units, where order parameter is scaled with $k_BT_c$, and 
    magnetic field is measured in tesla.

    Parameters
    ----------
    p : float, array-like
        Pressure in bar.

    Returns
    -------
    Value of $g_H$, evaluated with pressure-dependent $F_0^a$ and $T_c$.

    """
    return h3c.gH0/((1 + h3d.F0a_poly(p))*(Tc_mK(p)))**2

def gz(p):
    """
    Material parameter for PHA magnetic energy, linear in magnetic field H, in 
    dimensionless units, where order parameter is scaled with $k_BT_c$, and 
    magnetic field is measured in tesla.

    Parameters
    ----------
    p : float or arraay-like
        Pressure in bar.

    Returns
    -------
    Value of $g_z$, evaluated with pressure-dependent $F_0^a$ and $T_c$.

    """
    bn = beta_norm_asarray(1, p)
    bn5 = bn[4]
    bn24 = bn[1] + bn[3]
    return (-bn5/bn24)*h3c.lambda_A1/Tc_mK(p)

def H_scale(t, p):
    r"""
    Magnetic field scale implied by material parameters $\alpha$ and $g_H$, 
    through $H_s^2 \equiv \sqrt{-\alpha/g_H}$).

    Parameters
    ----------
    t : float
        reduced temperature.
    p : float or arraay-like
        Pressure in bar.

    Returns
    -------
    float or arraay-like (same shape as p)
        Magnetic field in tesla.

    """
    return np.sqrt(-alpha_norm(t)/gH(p))

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


