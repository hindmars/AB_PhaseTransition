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

def get_setting(name):
    """Gets value of a globally defined variable.
    """
    return globals()[name]

def report_setting(name):
    """Reports value of a globally defined variable.
    """
    xval = globals()[name]
    print("he3_tools:", type(xval), "variable " + name + " set to", xval)
    return

def set_default(name, xval):
    """Sets value of a globally defined variable.
    """
    set_name = name.replace("DEFAULT_", "SET_")
    if xval in globals()[set_name]:
        globals()[name] = xval
    else:
        raise ValueError("error: " + xval + " not in " + set_name)
    report_setting(name)
    return 

def b_wc(n):
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
        b = h3c.b_wc_list[n-1]
    else:
        raise ValueError("b_wc: n should be 1, 2, 3, 4 , or 5")
    return b

def b_wc_array():
    """
    Parameters
    ----------

    Returns
    -------
    b : float (or python default, usually is double)
        Normalised b as (5,) shape array, for weak coupling app.

    """
    b = h3c.b_wc_list
    return np.array(b)

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
    """Converts reduced temperature ($T/T_c$) to temperature in mK.
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
    Melting pressure , Greywall 1986 (equation A1). T in mK, p in bar.
    Uses Tian, Smith, Parpia 2022 PLTS to Greywall temperature converter.
    """
    if isinstance(T_mK, float) or isinstance(T_mK, int):
        T = T_mK
    else:
        T = T_mK.copy()
    if DEFAULT_T_SCALE == "PLTS":
        T = T_PLTStoG(T_mK)
    # return h3d.p_A_bar * np.ones_like(T_mK)
    # Greywall data down to 0.9mK onlyq
    T = np.maximum(T, 0.9)
    return h3d.Pmelt_poly_Greywall(T)/T**3 + h3d.p_A_bar

def T_melt_PLTS(p):
    """Returns melting temperataure (PLTS) as a function of pressure. 
    Uses Tian, Smith, Parpia NIMS 2022 interpolation polunomials.
    """
    if isinstance(p, float) or isinstance(p, int):
        pa = np.array([p])
    else:
        pa = p
        
    dp_mbar = (pa - h3d.p_A_bar)*1e3 

    T_melt_mK = np.ones_like(pa)*np.nan
    T_melt_mK[pa > 29.3] = h3d.Tmelt_poly_PLTS_hi(dp_mbar[pa > 29.3])
    T_melt_mK[pa > 35.2] = h3d.Tmelt_poly_PLTS_lo(dp_mbar[pa > 35.2])
    
    if isinstance(p, float) or isinstance(p, int):
        T_melt_mK = T_melt_mK[0]
    
    return T_melt_mK

def T_GtoPLTS(TG):  
    """Temperature scale convertor, Greywall to PLTS, ninth order polynomial.
    Between 0.9 mK and 100 mK, uses Tian, Smith, Parpia NIMS 2022 interpolation 
    polynomials. For T > 100 mK, assumes linear realtio consistent with TSP value 
    at 100 mK.    
    """
    if isinstance(TG, float) or isinstance(TG, int):
        TGa = np.array([TG])
    else:
        TGa = TG
        

    T_PLTS_mK = np.ones_like(TGa)*np.nan
    T_PLTS_mK[TGa <= 5.6] = h3d.GtoPLTS_low_poly(TGa[TGa <= 5.6])
    T_PLTS_mK[TGa > 5.6] = h3d.GtoPLTS_high_poly(TGa[TGa > 5.6])
    T_PLTS_mK[TGa > 100] = TGa[TGa > 100] - 100.0 + h3d.GtoPLTS_high_poly(100)
    
    if isinstance(TG, float) or isinstance(TG, int):
        T_PLTS_mK = T_PLTS_mK[0]

    return T_PLTS_mK

def T_PLTStoG(TPLTS):
    """Temperature scale converter, invert Greywall to PLTS, ninth order polynomial 
    Works between 0.9 mK and 100 mK.
    Uses Tian, Smith, Parpia NIMS 2022 interpolation polunomials.        
    """
    return h3d.PLTStoG9_high_poly(TPLTS)

def npart(p):
    """Particle density at pressure p bar ($nm^{-3}$).
    """
    # return np.interp(p, p_nodes, np_data)
    return np.interp(p, h3d.p_nodes, h3d.np_data)

def mstar_m(p):
    """Effective mass ratio at pressure p bar.
    """
    return np.interp(p, h3d.p_nodes, h3d.mstar_m_data)

def vf(p):
    """Fermi velocity at pressure p bar (m/s).
    """
    return np.interp(p, h3d.p_nodes, h3d.vf_data)

def pf(p):
    """Fermi momentum at pressure p bar (SI units, kg m/s)
    """
    return h3c.mhe3_kg * mstar_m(p) * vf(p)

def muf(p, units='SI'):
    """Chemical potential at pressure p bar (Fermi energy at zero temperature).
    units = 'SI' (default) SI units (J)
    units = 'eV' (electron-volts)
    
    """
    ef = 0.5*h3c.mhe3_kg * mstar_m(p) * vf(p)**2
    if units == 'eV':
        ef /= h3c.c.e
    return ef

def xi0(p):
    """Zero T Cooper pair correlation length at pressure p bar (nm).
    """
    return np.interp(p, h3d.p_nodes, h3d.xi0_data)

def F0a(p):
    """Landau parameter $F_0^a$."""
    return h3d.F0a_poly(p)

def xi(t, p):
    """Ginzburg Landau correlation length at pressure p bar (nm).
    """
    return h3c.xiGL_const*xi0(p)/(-alpha_norm(t))**0.5

def xi_delta(t, p):
    """BCS correlation length at pressure p bar (nm).
    """
    return h3c.xiGL_const*xi0(p)/(-alpha_bcs(t))**0.5
    
def N0(p):
    """Density of states at Fermi surface, units nm$^{-3}$ J$^{-1}$.
    """
    return npart(p) * 1.5 / (h3c.mhe3_kg * mstar_m(p) * vf(p)**2)

def Gi(p):
    """Ginzrurg number at pressure p.
    
    Gi = N(0) * xi0(p)**3 * kB * T_c(p)
    
    N(0) - density of states at Fermi surface
    xi0  - Cooper pair correlation length
    T_c  - normal/superfluid critical temperature
    
    """
    
    return N0(p) * xi0(p)**3 * h3c.kB * Tc_mK(p)*1e-3

def tGL(p):
    r"""
    Ginzburg-Landau time, defined as 
    $$
    t_{\rm GL} = \xi_{\rm GL}(p)/ v_F(p)
    $$
    where $\xi_{\rm GL}$ is the zero-temperature coherence length and $v_F$ is 
    the Fermi velocity.
    
    Parameters
    ----------
    p : float, int or numpy.ndarray
        Pressure in bar.

    Returns
    -------
    float, int or numpy.ndarray
        GL time.

    """
    return xi(0, p)/vf(p)

def tauN0(t, p):
    """
    Mean free timescale of quasiparticles (a.k.a. quasiparticle lifetime).
    Uses simple approximation for $\tau_N^0$ gicen in Vollhardt and Woelfle, p36,
    crediting Wheatley 1978. DEPRECATED in favour of tau_qp

    Parameters
    ----------
    t : float, np.ndarray
        Reduced temperature.
    p : float, np.ndarray
        Pressure in bar.
    
    Only one of t, p is allowed to be an array.

    Returns
    -------
    float, nd.array
        Quasiparticle lifetime in seconds.

    """
    
    return 0.3e-6 /T_mK(t, p)**2 

def tau_qp(t, p):
    """
    Mean free scattering time of quasiparticles (a.k.a. quasiparticle lifetime).
    Uses Greywall 1984 data, and enforces Greywall temparature.

    Units: s
    
    Parameters
    ----------
    t : float, np.ndarray
        Reduced temperature.
    p : float, np.ndarray
        Pressure in bar.
    
    Only one of t, p is allowed to be an array.

    Returns
    -------
    float, nd.array
        Quasiparticle lifetime in seconds.

    """
    p_data = h3d.data_Gre86_therm_cond[:,0]
    tauT2_data = h3d.data_Gre86_therm_cond[:,5]
    
    # Greywall is in sec/K^2, same as musec/mK^2
    tauT2_interp = np.interp(p, p_data, tauT2_data)
    print(p, tauT2_interp)
    
    if get_setting('DEFAULT_T_SCALE') != 'Greywall':
        tmp_set = get_setting('DEFAULT_T_SCALE')
        set_default('DEFAULT_T_SCALE', 'Greywall')
        T = T_mK(t, p)
        set_default('DEFAULT_T_SCALE', tmp_set)
    else:
        T = T_mK(t, p)
    print(p, tauT2_interp, T)
    
    return tauT2_interp /T**2 * 1e-6

def mfp0(t, p):
    """
    Mean free path of quasiparticles.
    Uses simple approximation for $\tau_N^0$ gicen in Vollhardt and Woelfle, p36,
    crediting Wheatley 1978.

    Parameters
    ----------
    t : float, np.ndarray
        Reduced temperature.
    p : float, np.ndarray
        Pressure in bar.
    
    Only one of t, p is allowed to be an array.

    Returns
    -------
    float, nd.array
        Quasiparticle mean free path in nm.

    """
    
    return vf(p)*tauN0(t, p)*1e9

def mfp_qp(t, p):
    """
    Mean free path of quasiparticles.  Uses Greywall data for $\tau_qp$. 
    
    Units: nm

    Parameters
    ----------
    t : float, np.ndarray
        Reduced temperature.
    p : float, np.ndarray
        Pressure in bar.
    
    Only one of t, p is allowed to be an array.

    Returns
    -------
    float, nd.array
        Quasiparticle mean free path in nm.

    """
    
    return vf(p)*tau_qp(t, p)*1e9

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

    """ Natural free energy density scale, units Joule per nm3 .
    """
    # return (1/3) * N0(p) * (2 * np.pi * kB * T_mK(1, p) * 1e-3)**2
    return (1/3) * N0(p) * (h3c.kB * Tc_mK(p) * 1e-3)**2
    
def delta_b(p, n):
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
        db = delta_b_interp(p, n)
    
    if DEFAULT_SC_ADJUST:
            db *= np.exp(-sc_adjust_fun(p))
    return db

def delta_b_asarray(p, n_list=range(1,6)):
    """
    Strong coupling corrections to material parameters, in units of 
    the modulus of the first weak coupling parameter, supplied as a (1,5) array.
    """
    delta_b_list = [ delta_b(p, n) for n in n_list]
    return np.array(delta_b_list)

def delta_b_interp(p, n): 
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

def delta_b_polyfit(p, n): 
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
    b = b_wc(n)
    # else:
    #     raise ValueError("beta_norm: n must be between 1 and 5")
    return h3c.beta_const*(b + t * delta_b(p, n))


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
    """Normalised free energy in a given phase.
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

def t_AB(p, low_pressure_nan=True):
    """ AB transition temperature at pressure p, normalised to Tc.
    """
    t_ab_val = (1/3)/ (delta_b(p, 1) 
                       + (delta_b(p, 3) 
                       - 2*delta_b(p, 4) 
                       - 2*delta_b(p, 5))/3) 
    
    if low_pressure_nan:
        if isinstance(t_ab_val, np.ndarray):
            t_ab_val[t_ab_val > 1] = np.nan
        elif t_ab_val > 1:
            t_ab_val = np.nan
            
    return  t_ab_val

def tAB(p, low_pressure_nan=True):
    """Synonym for t_AB.
    """
    return t_AB(p, low_pressure_nan)

def TAB_mK(p):
    """ AB transition temperature at pressure p, in mK
    """
    return tAB(p) * Tc_mK_expt(p)

# specific heats

def C_V_normal(t, p):
    """
    Normal phase specific heat in units J / K / nm$^3$, with strong coupling 
    corrections included.  Uses formula (2.13) from Vollhardt & Woelfle.
    $$
    C_V = \frac{\pi^2}{3} N_F k_B^2 T
    $$

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Normal phase specific heat at temperature $t$ and pressure p.

    """
    
    return (np.pi**2/3) * h.N0(p) * h.kB**2 *( h.Tc_mK(p)/1000) * t    

def C_V_normal(t, p):
    """
    Normal phase specific heat in units J / K / nm$^3$, with strong coupling 
    corrections included.  Uses formula (2.13) from Vollhardt & Woelfle.
    $$
    C_V = \frac{\pi^2}{3} N_F k_B^2 T
    $$

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Normal phase specific heat at temperature $t$ and pressure p.

    """
    
    # return (np.pi**2/3) * N0(p) * h3c.kB**2 *(Tc_mK(p)/1000) * t    
    return  np.pi**2 * f_scale(p) / (Tc_mK(p)/1000) * t


def delta_C_V_Tc(p, phase):
    """
    Jump in specific heat at superfluid phase transition in units J / K / nm$^3$, 
    with strong coupling corrections included.  Uses formula (3.78) from 
    Vollhardt & Woelfle.
    $$
    \Delta C_V = \frac{1} N(0) \frac{\partial}{\partial T} \langle \Delta_{\bf k}^\dagger \Delta_{\bf k} \rangle_{\hat{\bf k}}
    $$

    Parameters
    ----------
    p : float, int, numpy.ndarray
        Pressure in bar.
    phase : str
        Phase string, e.g. "A", "B".
        

    Returns
    -------
    float or numpy.ndarray
        Specific heat jump at pressure p.

    """
    return  0.25 * f_scale(p) / (Tc_mK(p)/1000 * beta_phase_norm(1, p, phase))

def C_V(t, p, phase):
    """
    Specific heat in given phase, with strong coupling 
    corrections included.  Uses formula (3.77) from Vollhardt & Woelfle.
    $$
    C_V = C_N + \Delta C_V
    $$

    Units:  units J / K / nm$^3$

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Specific heat at temperature $t$ and pressure p.

    """
    return C_V_normal(t,p) + delta_C_V_Tc(p, phase)

def kappa_0(t, p):
    """
    Gas-kinetic expression for thermal conductivity, V&W eqn 2.40.  
    
    Units: J / ns / nm / K

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Gas-kinetic thermal conductivity at temperature $t$ and pressure p.

    """
    # C_V is in J / K / nm^3, vf is in m/s (also nm/ns), tau_qp is in seconds
    return (1/3) * C_V_normal(t, p) * (vf(p))**2 * tau_qp(t, p) * 1e9

def kappa(t, p):
    """
    Experimental thermal conductivity, Greywall 1984  
    
    Units:    J / ns / nm / K

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Measured thermal conductivity at temperature $t$ and pressure p.

    """
    p_data = h3d.data_Gre86_therm_cond[:, 0]
    kappaT_data = h3d.data_Gre86_therm_cond[:, 4]
    
    # Greywall is in erg/sec/cm, convert to J/ns/nm 
    kappaT_interp = np.interp(p, p_data, kappaT_data) * 1e-7/(1e9 * 1e9 * 1e-2)

    if get_setting('DEFAULT_T_SCALE') != 'Greywall':
        tmp_set = get_setting('DEFAULT_T_SCALE')
        set_default('DEFAULT_T_SCALE', 'Greywall')
        T = T_mK(t, p)*1e-3
        set_default('DEFAULT_T_SCALE', tmp_set)
    else:
        T = T_mK(t, p)*1e-3

    print(p, T, kappaT_interp)

    return kappaT_interp/T

def gamma_c(p):
    """
    Order parameter damping rate at $T_c$. Units are ns$^{-1}$

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Order parameter damping rate at $T_c$

    """
    return (16 * h3c.kB * Tc_mK(p)*1e-3) / (np.pi * h3c.hbar) * 1e-9

def thermal_diffusivity(t, p):
    r"""
    Thermal diffusivity, defined as
    $$
    D_t = \kappa / C_V,
    $$
    where $\kappa$ is thermal conductivity and $C_V$ is volumetric heat capacity.
    
    Units: nm$^2$ / ns

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Thermal diffusivity.

    """
    return kappa(t, p)/C_V_normal(t, p)

def thermal_diffusivity_norm(t, p):
    r"""
    Dimensionless thermal diffusivity, defined as
    $$
    D_t t_{\rm GL} / \xi_{\rm GL}^2
    $$
    where $\xi_{\rm GL}$ is the zero-temperature GL coherence length and 
    $t_{\rm GL} = \xi_{\rm GL}/v_f$ is the GL time.    

    Parameters
    ----------
    t : float, int, numpy.ndarray
        Reduced temperature, $T/T_c$.
    p : float, int, numpy.ndarray
        Pressure in bar.

    Only one or other of t and p can be an array.

    Returns
    -------
    float or numpy.ndarray
        Thermal diffusivity.

    """
    return thermal_diffusivity(t, p) * tGL(p) / xi(0, p)**2

def diffusion_relaxation_length(p):
    r"""
    A length scale derived from the thermal diffusivity (dimesnions L$^2$ / T) 
    and the order parameter relaxation time at the critical temperature 
    $\gamma_c$, dimensions T${-1}$.

    Parameters
    ----------
    p : float, int, numpy.ndarray
        Pressure in bar.

    Returns
    -------
    float or numpy.ndarray
        Diffusion relaxation length.

    """
    return (thermal_diffusivity(1, p) / gamma_c(p)  )**0.5

def hot_blob_length_scale(p, phase='A'):
    """
    Length scale of a normal region produced by 1 eV energy, injected into 
    superfluid in given phase at zero temperature, defined as 
    $$
    L = (C_V T_c/ 1 {\rm eV})^{-1/3}
    $$
    Units: nm.
    
    Radius of hot blob would contain geometric factors from region shape, and 
    algebraic factors from integration of temperature between $0$ and $T_c$.

    Parameters
    ----------
    p : float, int, numpy.ndarray 
        Pressure in bar.

    Returns
    -------
    type(p)
        Length scale of the heated region, in nm..

    """
    C_V_Tc = C_V_normal(1, p) + delta_C_V_Tc(p, phase)
    return (Tc_mK(p)/1000 * C_V_Tc/h3c.c.e)**(-1/3)

def hot_blob_size(p, t, Q_eV, phase='A'):
    """
    Size of a normal region produced by Q_eV energy, injected into 
    superfluid in given phase at reduced temperature t. 
    
    Units: nm.
    
    Radius of hot blob would contain geometric factors from region shape, and 
    algebraic factors from integration of temperature between $0$ and $T_c$.

    Parameters
    ----------
    p : float, int, numpy.ndarray 
        Pressure in bar.

    Returns
    -------
    type(p)
        Length scale of the heated region, in nm..

    """
    return hot_blob_length_scale(p, phase='A') * (Q_eV/(1 - t))**(1/3)

def hot_blob_quench_time(p, t, Q_eV=1.0, phase='A'):
    """
    Quench time of a normal region produced by Q_eV energy, injected into 
    superfluid in given phase at reduced temperature t. 
    
    Units: ns.
    
    Parameters
    ----------
    p : float, int, numpy.ndarray 
        Pressure in bar.

    Returns
    -------
    type(p)
        Quench time for hot blob, in ns..

    """
    return hot_blob_size(p, t, Q_eV, phase='A')**2/thermal_diffusivity(t, p) *1/(1-t)

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


