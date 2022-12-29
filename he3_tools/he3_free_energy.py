#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 17:18:50 2022

Functions for returning bulk free energy and its derivatives.

@author: hindmars
"""

import numpy as np
import he3_props as h3p
import he3_bases as h3b
import he3_matrix as h3m

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


# def U(A, alpha_norm, beta_norm_arr):
def U(A, *args):
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
    al, bn, _, _ = args_parse(*args)

    dim = A.ndim
    
    A_T = np.swapaxes(A, dim-2, dim-1)
    A_C = np.conj(A)
    A_H = np.conj(A_T)

    Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
    Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
    Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
    Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
    Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
    Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
    return (Un0 + Un1 + Un2 + Un3 + Un4 + Un5).real


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
    al, bn, _, _ = args_parse(*args)

    dim = A.ndim
    
    A_T = np.swapaxes(A, dim-2, dim-1)
    A_C = np.conj(A)
    A_H = np.conj(A_T)

    Un0 = al * h3m.tr( np.matmul(A , A_H) )
    
    Un1 = bn[0] *  h3m.tr( np.matmul(A , A_T)) * h3m.tr(np.matmul(A_C , A_H) ) 
    Un2 = bn[1] *  h3m.tr( np.matmul(A , A_H) )**2
    Un3 = bn[2] *  h3m.tr( np.matmul(np.matmul(A , A_T) , np.matmul(A_C , A_H) )  )
    Un4 = bn[3] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A   , A_H) )  )
    Un5 = bn[4] *  h3m.tr( np.matmul(np.matmul(A , A_H) , np.matmul(A_C , A_T) )  )
    return np.array([Un0, Un1, Un2, Un3, Un4, Un5])


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

    Returns
    -------
    float or ndarray
        Normalised derivative of bulk free energy.

    """
    al, bn, _, _ = args_parse(*args)

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
    return dUn0 + dUn1 + dUn2 + dUn3 + dUn4 + dUn5


def line_section(X, D, t, p, scale=None, n=500):
    """
    Retuens a line section in order parameter space, along with the free energy 
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
  
    if isinstance(X, str):
        delta_X = h3p.delta_phase_norm(t, p, X)
        X = h3b.D_dict[X]
        if isinstance(D, str):
            D = h3b.D_dict[D] - X
            D = D/h3m.norm(D)

    A_XD = np.multiply.outer( np.ones_like(v) , X)*delta_X + np.multiply.outer( v , D)*scale
    U_XD = U(A_XD, (t-1), h3p.beta_norm_asarray(t, p) )
    
    return v, A_XD, U_XD