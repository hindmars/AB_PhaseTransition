#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 17:18:50 2022

Functions for returning bulk free energy and its derivatives, as function of 
order paraameter. Physical conditions can be specified either as (t,p,H), or a 
set of material paaraameters alpha, beta, gH, gz and H.  Magnetic field H can be 
either a float or a vector. If H is a float, the magnetic field is assuming in 
the direction he3_bases.e[0] = np,array([0,0,1]).

@author: hindmars
"""

import numpy as np
import he3_props as h3p
import he3_bases as h3b
import he3_matrix as h3m
import scipy.optimize as spo

def args_parse(*args):
    n_el = len(args)
    if n_el == 2:
        if isinstance(args[1], np.ndarray):
            t = None
            p = None
            al = args[0]
            bn = args[1]
        else:
            t, p = args
            al = h3p.alpha_norm(t)
            bn = h3p.beta_norm_asarray(t, p)
    else:
        t = None
        p = None
        al = args[0][0]
        bn = args[0][1:]
    return al, bn, t, p

def Uargs_parse(*args):
    """
    Parses arguments passed to free energy to determine material parameters.  
 
    Parameters
    ----------
    *args : floats or arraays of floats
       *(t,p) ; (float, float) - reduced temperature and pressure [bar].
       *(t,p,H) ; (float, float, float) - reduced temperature, pressure [bar], and 
               magnetic field [tesla], assumed pointing in z direction
       *(t,p,H_arr) ; (float, float, ndarray) - reduced temperature, pressure [bar], and 
               magnetic field vector [tesla] with shape (3,).
       *(alpha, beta_arr) ; (float, ndarray) - material parameters alpha and beta_a. 
               beta_arr must have shape (5,).
       *(alpha, beta_arr, gH, gz, H) ; (float, ndarray, 3*float) - material parameters 
               alpha, beta_a, gH, gz and magnetic field [tesla] H, assumed || to z.
               beta_arr must have shape (5,).
       *(alpha, beta_arr, gH, gz, H_arr) ; (float, ndarray, 2*float, ndarray) - material parameters 
               alpha, beta_a, gH, gz and magnetic field [tesla] H, shape (3,).
               beta_arr must have shape (5,1).

    Raises
    ------
    ValueError
        If number of args is not 2, 3 or 5.
        If second arg is array, and number of args is not 2 or 5.

    Returns
    -------
    al : float
        Material parameter alpha.
    bn : ndarray, shape (5,1).
        Material parameters beta.
    t : float
        Reduced temperature, $t = T/T_c(p)$.
    p : float
        Pressure [bar].

    """
    n_el = len(args)
    gH = 0.
    gz = 0.
    H = 0.
    if n_el not in [2,3,5]:
        raise ValueError('Uargs_parse: error: number of args should be 2, 3 or 5')
    if isinstance(args[1], np.ndarray):
        t = None
        p = None
        al = args[0]
        bn = args[1]
        gH = 0
        gz = 0
        H = 0*h3b.e[0]
        if n_el == 5:
            gH = args[2]
            gz = args[3]
            H = args[4]
        elif n_el != 2:
            raise ValueError('Uargs_parse: second argument is an array, expecting 2 or 5 arguments')
    else:
        t, p = args[0:2]
        al = h3p.alpha_norm(t)
        bn = h3p.beta_norm_asarray(t, p)
        if n_el == 3:
            gH = h3p.gH(p)
            gz = h3p.gz(p)
            H = args[2]
        
    if isinstance(H, float) or isinstance(H, int):
        H = H*h3b.e[0]

    # else:
    #     t = None
    #     p = None
    #     al = args[0][0]
    #     bn = args[0][1:]
    return al, bn, gH, gz, t, p, H

# def U(A, alpha_norm, beta_norm_arr):
def U(A, *args):
    """
    Bulk free energy for superfluid He3.  See Uargs_parse for allowed combinations 
    of parameters.

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter.
    alpha_norm : float
        alpha parameter, normalised.
    beta_norm_arr : ndarray, shape (5,1)
        beta parameters, normalised.
    gH_norm : float, optional
        gH parameter, normalised.
    gz_norm : float, normalised
        gH parameter, normalised.
    t : float, optional
        reduced temperature
    p : float, optional
        pressure [bar]
    H_T : float or ndarray shape (3,), optional
        H field, (in tesla )

    Returns
    -------
    float or ndarray
        Normalised bulk free energy.

    """
    # al, bn, _, _ = args_parse(*args)
    al, bn, gH, gz, t, p, H = Uargs_parse(*args)
    dim = A.ndim
    A_T = np.swapaxes(A, dim-2, dim-1)
    # A_C = np.conj(A)
    A_H = np.conj(A_T)

    AA_H = np.matmul(A , A_H)
    AA_T = np.matmul(A , A_T)
    
    # Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
    # Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
    # Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
    # Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
    # Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
    # Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
    
    Un0 = al * h3m.tr( AA_H )
    
    Un1 = bn[0] *  h3m.tr( AA_T) * np.conj(h3m.tr( AA_T))
    Un2 = bn[1] *  h3m.tr( AA_H )**2
    Un3 = bn[2] *  h3m.tr( np.matmul( AA_T , np.conj(AA_T))  )
    Un4 = bn[3] *  h3m.tr( np.matmul( AA_H , AA_H )  )
    Un5 = bn[4] *  h3m.tr( np.matmul( AA_H , np.conj(AA_H) )  )
    fHn2 = gH * np.dot(np.matmul(AA_H, H), H)
    fHn1 = -1j * gz * np.dot(h3m.levi_civita3_matrix2vector(AA_H), H)
    # print((Un0 + Un1 + Un2 + Un3 + Un4 + Un5 + fHn2 + fHn1).imag)
    return (Un0 + Un1 + Un2 + Un3 + Un4 + Un5 + fHn2 + fHn1).real
    # return (Un0 + Un1 + Un2 + Un3 + Un4 + Un5).real


# def U(A, alpha_norm, beta_norm_arr):
def U_terms(A, *args):
    """
    Bulk free energy for superfluid He3

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
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
    # al, bn, _, _ = args_parse(*args)
    al, bn, gH, gz, t, p, H = Uargs_parse(*args)

    dim = A.ndim
    
    # A_T = np.swapaxes(A, dim-2, dim-1)
    # A_C = np.conj(A)
    # A_H = np.conj(A_T)

    # Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
    # Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
    # Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
    # Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
    # Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
    # Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )

    A_T = np.swapaxes(A, dim-2, dim-1)
    A_H = np.conj(A_T)

    AA_H = np.matmul(A , A_H)
    AA_T = np.matmul(A , A_T)
    
    Un0 = al * h3m.tr( AA_H )
    
    Un1 = bn[0] *  h3m.tr( AA_T) * np.conj(h3m.tr( AA_T))
    Un2 = bn[1] *  h3m.tr( AA_H )**2
    Un3 = bn[2] *  h3m.tr( np.matmul( AA_T , np.conj(AA_T))  )
    Un4 = bn[3] *  h3m.tr( np.matmul( AA_H , AA_H )  )
    Un5 = bn[4] *  h3m.tr( np.matmul( AA_H , np.conj(AA_H) )  )
    fHn2 = gH * np.dot(np.matmul(AA_H, H), H)
    fHn1 = - 1j * gz * np.dot(h3m.levi_civita3_matrix2vector(AA_H), H)    
    return np.array([Un0, Un1, Un2, Un3, Un4, Un5, fHn2, fHn1])


def dU_dA(A, *args):
    """
    Derivative of bulk free energy for superfluid He3.')

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...3,3)
        order parameter.
    alpha_norm : float
        alpha material parameter, dimensionless.
    beta_norm_arr : ndarray, shape (5,1)
        beta material parameters, dimensionless.
    gH_norm : float, optional
        gH parameter, normalised.
    gz_norm : float, normalised
        gH parameter, normalised.
    t : float, optional
        reduced temperature
    p : float, optional
        pressure [bar]
    H_T : float or ndarray shape (3,), optional
        H field, (in tesla )

    Returns
    -------
    float or ndarray
        Normalised derivative of bulk free energy.

    """
    # al, bn, _, _ = args_parse(*args)
    al, bn, gH, gz, t, p, H = Uargs_parse(*args)

    dim = A.ndim
    
    A_T = np.swapaxes(A, dim-2, dim-1)
    A_C = np.conj(A)
    A_H = np.conj(A_T)

    dUn0 = al * A
        
    dUn1 = 2 * bn[0] *  np.matmul(np.multiply.outer(h3m.tr( np.matmul(A , A_T)), h3b.id3 ) , A_C)
    dUn2 = 2 * bn[1] *  np.matmul(np.multiply.outer(h3m.tr( np.matmul(A , A_H)), h3b.id3 ) , A )
    dUn3 = 2 * bn[2] *  np.matmul(A , np.matmul(A_T , A_C))
    dUn4 = 2 * bn[3] *  np.matmul(A , np.matmul(A_H   , A))
    dUn5 = 2 * bn[4] *  np.matmul(A_C , np.matmul(A_T , A))
    
    dfHn2 = gH * np.multiply.outer(np.matmul(H, A), H)
    dfHn1 = -1j*gz * np.swapaxes(np.matmul(A_T, h3m.levi_civita3_vector2matrix(H)), dim-2, dim-1)
    # print(dfHn2)
    # print(dfHn1)
    return dUn0 + dUn1 + dUn2 + dUn3 + dUn4 + dUn5 + dfHn2 + dfHn1


def line_section(X, D, t, p, path='linear', scale=None, norm_preserve=False, n=500):
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
        delta_X = h3p.delta_phase_norm(t, p, X)
        X = h3b.D_dict[X]
        A_0 = X * delta_X
        if isinstance(D, str):
            delta_Y = h3p.delta_phase_norm(t, p, D)
            Y = h3b.D_dict[D]
            A_1 = Y * delta_Y
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
    U_XD = U(A_XD, t, p)
    
    return v, A_XD, U_XD

def get_A_ext(A_init, *args):
    """
    

    Parameters
    ----------
    A_init : complex numpy array shape (3,3)
        Order parameter initial guess.
    *args : tuple
        See Uargs_parse.
        
    Find extremum of potential U using Newton-Krylov method.

    Returns
    -------
    complex numpy array shape (3,3)
        Local extremum of potential nearest A_init.

    """
    

    def F(A):
        return dU_dA(A, *args)
    
    try:
        A = spo.newton_krylov(F, A_init)
    except spo.nonlin.NoConvergence as e:
        A = e.args[0]
        print('get_U_min: No Convergence')
    
    return A
