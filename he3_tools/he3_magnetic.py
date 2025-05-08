#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Created on Sun Jan  1 12:27:28 2023

Functions for finding gaps and free energies in presnece of magnetic field. Gaps 
are for $\Delta_{\uparrow\uparrow}$. 

@author: hindmars
"""
import numpy as np
import he3_constants as h3c
import he3_data as h3d
import he3_props as h3p
import he3_bases as h3b
import he3_matrix as h3m
import he3_free_energy as h3fe
from scipy.optimize import newton_krylov
from scipy.optimize import bisect
# import scipy.optimize.nonlin
from scipy.optimize import NoConvergence

def uB_mag_norm(delta_ss, t, p, H):
    if isinstance(p, int):
        p = float(p)

    al = h3p.alpha_norm(t)
    bn = h3p.beta_norm_asarray(t, p)    
    gH = h3p.gH(p)
    gz = h3p.gz(p)

    x, y, z = delta_ss

    uB2 = (1/9)*(bn[0]*(2*x*y+z**2)**2 + bn[1]*(x**2+y**2+z**2)**2 +\
        (bn[2]+bn[4])*(2*x**2*y**2+z**4) + bn[3]*(x**4+y**4+z**4)) +\
            (1/3)*(al*(x**2+y**2+z**2) + gH*H**2*z**2 - gz*H*(x**2-y**2))
    return uB2

def Hc_T(t,p):
    """Critical magnetic field for vanishing Delta_up_down in B phases, using 
    GL free energy.  If DEFAULT_ALPHA is set to BCS, it uses the he3_tools fit 
    to the weak coupling B-phase gap parameter to compute $\alpha$.
    """
    al = h3p.alpha_norm(t)
    bn = h3p.beta_norm_asarray(t, p)    
    bn12 = bn[0] + bn[1]
    bn345 = bn[2] + bn[3] + bn[4]
    c = h3p.gH(p)
    return  np.sqrt(-al*bn345/(c*(2*bn12 + bn345)))

def tcB2(p, H):
    """Critical reduced temperature for vanishing Delta_up_down."""
    c = h3p.gH(p)
    bn = h3p.beta_norm_asarray(1, p)
    bn12 = bn[0] + bn[1]
    bn345 = bn[2] + bn[3] + bn[4]

    return 1 - (1 + 2*bn12/bn345)*c*H**2

def delta_B_mag_norm_approx(t, p, H):
    """ Gap vector for B2 phase. 
    (Approximate, linear expansion in particle-hole asymmetry parameter $g_z$ or $\eta$).
    Takes floats t, p, H and returns (3,) array, comprising $\Delta_{s_1s_2}$, 
    with $s_1,s_2$ = (up,up), (down, down) and (up,down).  
    In the B-phase, with the magnetic field along the z axis, these 
    are the diagonal elements of the order parameter matrix in a suitable basis.
    """
    if isinstance(p, int):
        p = float(p)
    if isinstance(t, int):
        t = float(t)

    al = h3p.alpha_norm(t)
    bn = h3p.beta_norm_asarray(t, p)    
    c = h3p.gH(p)
    eta = h3p.gz(p)
    
    bn2 = bn[1]
    bn4 = bn[3]
    bn12 = bn[0] + bn[1]
    bn345 = bn[2] + bn[3] + bn[4]
    bn12345 = bn12 + bn345
    
    if isinstance(p, float):

        if isinstance(t, float):

            if t < 1:
                if t > tcB2(p, H):
                    z = 0
                    xp = np.sqrt(3/2)*np.sqrt(-al/(2*bn12 + bn345))
                    # if xp != 0.0:
                    xm = 3*H*eta*xp/(6*al + 8*bn2*xp**2 + 8*bn4*xp**2)
                    # else:
                    # xm = 0.0
                else:
                    xp = np.sqrt(3/2)*np.sqrt((H**2*bn12*c - al*bn345)/
                                              (bn345*(3*bn12 + bn345)))
                    xm = 3*H*bn12345*eta*xp/(2*(3*H**2*bn2*c - 3*al*bn12345 + 3*al*bn2)) 
                    z =  np.sqrt(3/2)*np.sqrt((-c*H**2*(2*bn12 + bn345) - al*bn345)/
                                              (bn345*(3*bn12 + bn345)))
            else:
                xp, xm, z = (0, 0, 0)

        elif isinstance(t, np.ndarray):
            
            xp = np.zeros_like(t, dtype=float)
            xm = np.zeros_like(t, dtype=float)
            z = np.zeros_like(t, dtype=float)
    
            B1 = np.logical_and(t >= tcB2(p, H), t< 1.)
            B2 = t < tcB2(p, H)
            print(B1.shape, B2.shape)
            xp[B2] = np.sqrt(3/2)*np.sqrt((H**2*bn12[B2]*c - al[B2]*bn345[B2])/
                                      (bn345[B2]*(3*bn12[B2] + bn345[B2])))
            # xm[B2] = 3*H*bn12345[B2]*eta[B2]*xp[B2]/(2*(3*H**2*bn2[B2]*c[B2] - 3*al*bn12345[B2] + 3*al*bn2[B2])) 
            xm[B2] = 3*H*bn12345[B2]*eta*xp[B2]/(-6*H**2*bn2[B2]*c + 6*al[B2]*bn12345[B2] - 6*al[B2]*bn2[B2] 
                                                     + 8*bn12345[B2]*bn4[B2]*xp[B2]**2 + 8*bn2[B2]*bn345[B2]*xp[B2]**2)
            z[B2] =  np.sqrt(3/2)*np.sqrt((-c*H**2*(2*bn12[B2] + bn345[B2]) - al[B2]*bn345[B2])/
                                      (bn345[B2]*(3*bn12[B2] + bn345[B2])))
            print(xp[B1].shape, al.shape, bn12[B1].shape, bn345[B1].shape)
            xp[B1] = np.sqrt(3/2)*np.sqrt(-al[B1]/(2*bn12[B1] +bn345[B1]))
            xm[B1] = 3*H*eta*xp[B1]/(6*al[B1] + 8*bn2[B1]*xp[B1]**2 + 8*bn4[B1]*xp[B1]**2)
            
            
        else:
            raise ValueError('delta_B2_norm_approx: t must be float or ndarray.\n' + 
                             'Was {}.'.format(type(t)))
  


    elif isinstance(p, np.ndarray):
        xp = np.zeros_like(p, dtype=float)
        xm = np.zeros_like(p, dtype=float)
        z = np.zeros_like(p, dtype=float)

        B1 = np.logical_and(t >= tcB2(p, H), t< 1.)
        B2 = t < tcB2(p, H)
        xp[B2] = np.sqrt(3/2)*np.sqrt((H**2*bn12[B2]*c[B2] - al*bn345[B2])/
                                  (bn345[B2]*(3*bn12[B2] + bn345[B2])))
        # xm[B2] = 3*H*bn12345[B2]*eta[B2]*xp[B2]/(2*(3*H**2*bn2[B2]*c[B2] - 3*al*bn12345[B2] + 3*al*bn2[B2])) 
        xm[B2] = 3*H*bn12345[B2]*eta[B2]*xp[B2]/(-6*H**2*bn2[B2]*c[B2] + 6*al*bn12345[B2] - 6*al*bn2[B2] 
                                                 + 8*bn12345[B2]*bn4[B2]*xp[B2]**2 + 8*bn2[B2]*bn345[B2]*xp[B2]**2)
        z[B2] =  np.sqrt(3/2)*np.sqrt((-c[B2]*H**2*(2*bn12[B2] + bn345[B2]) - al*bn345[B2])/
                                  (bn345[B2]*(3*bn12[B2] + bn345[B2])))
        xp[B1] = np.sqrt(3/2)*np.sqrt(-al/(2*bn12[B1] +bn345[B1]))
        xm[B1] = 3*H*eta[B1]*xp[B1]/(6*al + 8*bn2[B1]*xp[B1]**2 + 8*bn4[B1]*xp[B1]**2)

    else:
        raise ValueError('delta_B2_norm_approx: p must be float or ndarray.\n' + 
                         'Was {}.'.format(type(p)))

        
    return np.array([xp+xm, xp-xm, z])
    
def f_B_mag_norm_approx(t, p, H, order=2):
    """ Bulk free energy for B2 phase. 
    (Approximate, gap vector evalulated in linear expansion in 
    particle-hole asymmetry parameter $g_z$ or $\eta$).
    Takes floats t, p, H and returns float.
    """
    
    if order == 0 or order==1:
        al = h3p.alpha_norm(t)
        bn = h3p.beta_norm_asarray(t, p)    
        c = h3p.gH(p)
        
        bn12 = bn[0] + bn[1]
        bn345 = bn[2] + bn[3] + bn[4]

        Hs = np.sqrt(- al/c)
        h = H/Hs
        
        f_B2 = al**2 * (-3*bn345 + 2*bn345*h**2 - (2*bn12+bn345)*h**4)/\
            (4*bn345*(3*bn12 + bn345))
        
    elif order == 2:
        f_B2 = uB_mag_norm(delta_B_mag_norm_approx(t, p, H), t, p, H)
        
    return f_B2

#
#
# A phases
#
#

def uA_mag_norm(delta_ss, t, p, H):
    if isinstance(p, int):
        p = float(p)

    al = h3p.alpha_norm(t)
    bn = h3p.beta_norm_asarray(t, p)    
    gz = h3p.gz(p)

    x, y, z = delta_ss

    u_mag = (1/4)*(bn[1]*(x**2+y**2)**2 + bn[3]*(x**2+y**2)**2 + 4*bn[4]*x**2*y**2)\
            + (1/2)*(al*(x**2+y**2+z**2) - gz*H*(x**2-y**2))
    return u_mag

def delta_A_mag_norm_approx(t, p, H):
    """ Gap vector for A1 and A2 phase. 
    (Approximate, linear expansion in particle-hole asymmetry parameter $g_z$ or $\eta$).
    t: floats, p: float or ndarray, H float and returns (3,)+p.shape array, 
    comprising $\Delta_{s_1s_2}$, with $s_1,s_2$ = (up,up), (down, down) and (up,down).  
    In the A-phase, with the magnetic field along the z axis, these 
    are the diagonal elements of the order parameter matrix in the standard basis.
    """
    if isinstance(p, int):
        p = float(p)
    al = h3p.alpha_norm(t)
    bn = h3p.beta_norm_asarray(t, p)
    eta = h3p.gz(p)
    
    bn5 = bn[4]
    bn24 = bn[1] + bn[3]
    bn245 = bn[1] + bn[3] + bn[4]
    
    if isinstance(p, float):
        
        x,y,z = (0.,0.,0.)
        
        if t < tcA1(p, H):
            if t < tcA2(p, H):
                x = np.sqrt(1/2)*np.sqrt(-al/bn245 - eta*H/bn5)
                y = np.sqrt(1/2)*np.sqrt(-al/bn245 + eta*H/bn5)
            else:
                x = np.sqrt(-(al-eta*H)/bn24)
                y = 0.
    elif isinstance(p, np.ndarray):
        x = np.zeros_like(p, dtype=float)
        y = np.zeros_like(p, dtype=float)
        z = np.zeros_like(p, dtype=float)
        
        A1 = np.logical_and(t < tcA1(p, H), t >= tcA2(p, H))
        A2 = t < tcA2(p, H)
        x[A2] = np.sqrt(1/2)*np.sqrt(-al/bn245[A2] - eta[A2]*H/bn5[A2])
        y[A2] = np.sqrt(1/2)*np.sqrt(-al/bn245[A2] + eta[A2]*H/bn5[A2])
        x[A1] = np.sqrt(-(al-eta[A1]*H)/bn24[A1])
    else:
        raise ValueError('delta_A1A2_norm_approx: p must be float or ndarray.\n' + 
                         'Was {}.'.format(type(p)))
        
    return np.array([x, y, z])
    
def f_A_mag_norm_approx(t, p, H):
    """ Bulk free energy for A1 A2 phases. 
    (Approximate, gap vector evalulated in linear expansion in 
    particle-hole asymmetry parameter $g_z$ or $\eta$).
    Takes floats t reduced temoerature), p (pressure in bar)), H (B-field in tesla) 
    and returns float (dimensionless free energy, units set by he3_tools.f_scale(p), 
    in J/m$^3$.
    """

    return uA_mag_norm(delta_A_mag_norm_approx(t, p, H), t, p, H)

def tcA1(p, H):
    """Critical reduced temoerature for appearance of A1 phase.
    Takes floats t reduced temoerature), p (pressure in bar)).
    """
    eta = h3p.gz(p)
    return 1 + eta*H

def tcA2(p, H):
    """Critical reduced temoerature for appearance of A2 phase.
    Takes floats t reduced temoerature), p (pressure in bar)).
    """
    eta = h3p.gz(p)
    bn = h3p.beta_norm_asarray(1, p)
    bn5 = bn[4]
    bn245 = bn[1] + bn[3] + bn[4]

    return 1 + (bn245/bn5)*eta*H


#
# Magnetic field at equal free energy for A2 and B2 phases
#   zero PHA approximation
#

def HAB_T(t, p):
    """Magnetic field in tesla at which free energies of A and B phases are equal.
    
    t : float
        reduced temperature
    p : float or ndarray
        pressure in bar
        
    Returns
    
    float or ndarray, same shape as p. (Magnetic field in tesla)
    
    """
    
    close_to_zero = -1.0e-15
    
    if isinstance(p, int):
        p = float(p)
        if t > 1:
            return np.nan

    if t > 1:
        return np.nan*np.ones_like(p)
    
    bn = h3p.beta_norm_asarray(t, p)
    
    bn245 = bn[1] + bn[3] + bn[4]
    bn12 = bn[0] + bn[1]
    bn345 = bn[2] + bn[3] + bn[4]
    
    disc = bn245*bn345*(6*bn12**2 - 6*bn12*bn245 + 5*bn12*bn345 - 2*bn245*bn345 + bn345**2)

    if isinstance(p, float):
        if disc >= 0:
            # hAB = np.sqrt((1/(3*L)) * (-1 - np.sqrt(disc)))
            h2 = (bn245*bn345 - np.sqrt(disc))/(bn245*(2*bn12 + bn345))
            # print('h2', h2)
            if h2 >= close_to_zero:
                hAB = np.sqrt(h2 - close_to_zero)
                # print('hAB', hAB)
            else:
                hAB = np.nan
        else:
            hAB = np.nan
    elif isinstance(p, np.ndarray):

        hAB = np.nan*np.ones_like(p)
        h2  = np.nan*np.ones_like(p)
        
        real = disc >= 0
        h2[real] = (bn245[real]*bn345[real] - np.sqrt(disc[real]))/\
            (bn245[real]*(2*bn12[real] + bn345[real]))
        h2pos = h2 >= 0
        hAB[h2pos] = np.sqrt(h2[h2pos])
    else:
        raise ValueError('HAB_T: p must be float or ndarray')        
            
    return hAB * h3p.H_scale(t,p)

def HAB0_T_expt(p):
    """Critical magnetic field for AB transition at $T=0$ according to fit given in 
    Inseob Hahn, Y.H. Tang, H.M. Bozler, C.M. Gould, Phuysica B, 194-196B, 815 (1994).
    """
    pA = 34.338 # bar
    B0 = 0.339  # tesla
    y = p/pA
    return  B0*(1 + 6.34*y - 2.50*y**2)/ (1 + 2.10*y)


def tAB_mag(p, H, t_init=None):
    """Reduced temperature of AB transition in magnetic field H/Tesla. 
    At low pressures will converge only for t_init very close to 1.0.
    """
    def F(t):
        return HAB_T(t, p) - H

    # if t_init is None:
    t_init = h3p.tAB(p, low_pressure_nan=False)

    print('tAB initial guess', t_init)

    try:
        # tAB = newton_krylov(F, t_init)
        tAB = bisect(F, 1e-5, min(1.0, t_init))
    except NoConvergence as e:
        tAB = np.nan
        print('tAB_mag: No Convergence', e.args[0])

    return tAB

def TAB_mK_mag(p, H):
    
    return tAB_mag(p,H)*h3p.Tc_mK(p)

def find_mag_phase(phase, *args):
    """
    Find order parameter at minimum of bulk free energy (A and B phases only).
    
    """
    al, bn, gH, gz, t, p, H_vec = h3fe.Uargs_parse(*args)
    H = np.sqrt(np.sum(H_vec**2))
    
    def F(A):
        return h3fe.dU_dA(A, *args)
    
    if phase == 'A':
        del_ss = delta_A_mag_norm_approx(t, p, H)
        A_init = np.outer((del_ss[0] * h3b.e[1] + del_ss[1]*h3b.e[-1])/np.sqrt(2), h3b.e[1])
    elif phase == 'B':
        del_ss = delta_B_mag_norm_approx(t, p, H)
        A_init = (del_ss[0]*np.outer(h3b.e[1],h3b.e[-1]) + \
                del_ss[1]*np.outer(h3b.e[-1],h3b.e[1]) + \
                    del_ss[2]*np.outer(h3b.e[0],h3b.e[0]))/np.sqrt(3)
    elif phase == 'planar':
        A_init = h3p.delta_phase_norm(t, p, phase)*h3b.D_dict[phase]
    else:
        raise ValueError('find_mag_phase: only phase = A or B implemented.')
    try:
        A = newton_krylov(F, A_init)
    except NoConvergence as e:
        A = e.args[0]
        print('find_mag_phase: no Convergence')
        
    return A

def f_phase_mag_norm(phase, *args):
    """
    Find value of  minimum of bulk free energy (A and B phases only).
    """

    al, bn, gH, gz, t, p, H = h3fe.Uargs_parse(*args)
    
    A_ext = find_mag_phase(phase, *args)
        
    return h3fe.U(A_ext, *args)

def line_section(X, D, t, p, H, path='linear', scale=None, norm_preserve=False, n=500):
    """
    Returns a line section in order parameter space, along with the free energy 
    along it, and the parameter of the line.

    Parameters
    ----------
    X : string or Complex array shape (3,3)
        Start point of line. If a string, specifies the inert phase.
    D : string or Complex array shape (3,3)
        If a string, the end point inert phase. If array, the direction of the line
    t : float
        Reduced temperature in terms of Tc(p).
    p : float
        Pressure in bar.
    H : float
        Magnetic field in tesla.
    scale : float, optional
        Multiplies direction matrix. The default is None, in which case the 
        scale is the weak coupling gap scale delta_wc(t).
    n : integer, optional
        Number of points on the line. The default is 500.

    Returns
    -------
    v : numpy.ndarray shape (n,)
        Line parameter, between 0 and v_max=1.5.
    A_XD : Complex array shape (n,3,3)
        Order parameter along line.
    U_XD : float array shape (n,)
        Free energy along line.

    """
    if scale is None:
        scale = h3p.delta_wc(t)
    v_max = 1.5
    v = np.linspace(0,1,n)*v_max
    ind1 = np.argmin(v-1.)
  
    if isinstance(X, str):
        # delta_X = find_mag_phase(X, t, p, H)
        # X = h3b.D_dict[X]
        # A_0 = X * delta_X
        A_0 = find_mag_phase(X, t, p, H)
        if isinstance(D, str):
            # delta_Y = find_mag_phase(D, t, p, H)
            # Y = h3b.D_dict[D]
            # A_1 = Y * delta_Y
            A_1 = find_mag_phase(D, t, p, H)
            # A_XD = np.multiply.outer( 1-v , X)*delta_X + np.multiply.outer( v , Y)*delta_Y
        else:
            A_1 = A_0 + D * scale
            # A_XD = np.multiply.outer( np.ones_like(v) , X)*delta_X + np.multiply.outer( v , D)*scale
    else:
        # diff = D - X
        # A_XD = np.multiply.outer( np.ones_like(v) , X) + np.multiply.outer( v , diff)
        A_0 = X
        A_1 = D

    if path == 'linear':
        A_XD = np.multiply.outer( 1-v , A_0) + np.multiply.outer( v , A_1)
    elif path == 'OC77':
        prod = A_0 * A_1
        zer0 = np.zeros_like(prod)
        diff = np.float64(np.equal( prod, zer0))
        common = np.float64(np.not_equal( prod, zer0))
                
        theta = v * np.pi/2
        
        # delta_interp = np.cos(theta)**2 * h.delta_A_norm(t, p) + np.sin(theta)**2 * h.delta_planar_norm(t, p)
        
        A_XD = np.multiply.outer(np.cos(theta)**2, A_0*common) \
            + np.multiply.outer(np.sin(theta)**2, A_1*common)  \
            + np.multiply.outer(np.cos(theta), A_0*diff) \
            + np.multiply.outer(np.sin(theta), A_1*diff)
    else:
        raise ValueError('line_section: path must be "linear" or "OC77"')

    if norm_preserve:
        target_norm = h3m.norm(A_XD[0])*(1-v) + h3m.norm(A_XD[ind1])*v
        A_XD = np.matmul(np.multiply.outer(target_norm/h3m.norm(A_XD), h3b.id3),  A_XD)
    
    # U_XD = U(A_XD, (t-1), h3p.beta_norm_asarray(t, p) )
    U_XD = h3fe.U(A_XD, t, p, H)
    
    return v, A_XD, U_XD

def surface_energy_AB_approx(t, p, H=0, sigma=0.95):
    """Surface energy of AB interface, uses sigma_AB = sigma * f_B * xi_GL
    """
    return sigma*np.abs(f_phase_mag_norm('B', t, p, H))*h3p.xi(t,p)

def critical_radius(t_in, p, H=0, sigma=0.95, dim=3):
    """Radius of thin-wall critical bubble, in nm, with magnetic field. 
    Units: nm

    Ideally will optionally use function to get surface tension. 
    Uses approximations surface energy."""
    # if isinstance(sigma_fun, float):
    # sigma_AB = sigma*np.abs(f_B_mag_norm_approx(t,p,H))*h3p.xi(t,p)
    
    try:
        t_iterator = iter(t_in)
        f_B_mag_norm = np.zeros_like(t_in)
        f_A_mag_norm = np.zeros_like(t_in)
        for n, t in enumerate(t_in):
            f_B_mag_norm[n] = f_phase_mag_norm('B', t, p, H)
            f_A_mag_norm[n] = f_phase_mag_norm('A', t, p, H)
            
    except TypeError:    
        t = t_in
        f_B_mag_norm = f_phase_mag_norm('B', t, p, H)
        f_A_mag_norm = f_phase_mag_norm('A', t, p, H)
    
    sigma_AB = sigma*np.abs(f_B_mag_norm)*h3p.xi(t_in, p)
    # elif isinstance(sigma_fun, np.ndarray):
        # sigma_AB = sigma_fun*np.abs(f_B_norm(t,p))*xi(t,p)
    
    return (dim-1)*sigma_AB/(f_A_mag_norm - f_B_mag_norm)

def critical_energy(t, p, H=0, sigma=0.95, dim=3):
    """Energy of thin-wall critical bubble, in nm, with magnetic field.
    Ideally will optionally use function to get surface tension. 
    Uses approximations for free energy and surface energy.
    
    Units: k_B T_c
    """
    f_B_mag_norm = f_phase_mag_norm('B', t, p, H)
    f_A_mag_norm = f_phase_mag_norm('A', t, p, H)
    
    r_c = critical_radius(t, p, H, sigma, dim)
    sigma_AB = surface_energy_AB_approx(t, p, H, sigma)
    
    Ec = 4*np.pi*(r_c**2*sigma_AB + (1/3)*r_c**3*(f_B_mag_norm - f_A_mag_norm) )
    
    return Ec * h3p.f_scale(p)/(h3c.kB * h3p.Tc_mK(p)*1e-3)

def delta_tAB_mag_expt(p, H):
    
    p_data = h3d.data_Tan91_mag_supp_B[:,0]
    DT_mK_DB2_data = h3d.data_Tan91_mag_supp_B[:,1]

    DT_mK = np.interp(p, p_data, DT_mK_DB2_data*H**2)

    # print(p_data,  DT_mK_DB2_data, DT_mK)
    # print(p,  DT_mK.shape, h3p.Tc_mK(p).shape)
    # print(h3p.Tc_mK(p))
    
    return DT_mK/(h3p.Tc_mK(p)*tcB2(p,H))
    
    