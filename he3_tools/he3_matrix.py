#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 17:05:58 2022

Functions for matrix operations on order parameter arrays. The matrix indices are 
the last two in the array.

@author: hindmars
"""

import numpy as np
import scipy.linalg as sl


def tr(A):
    """
    Trace of order parameter array on last two indices

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d = A.ndim    
    return np.trace(A, axis1=d-2,  axis2=d-1 )


def trans(A):
    """
    Transpose of order parameter array on last two indices

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d = A.ndim    
    return np.swapaxes(A, d-2,  d-1 )


def hconj(A):
    """
    Hermitian conjugate of order parameter array on last two indices

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return trans(A).conj()


def mult_hconj(A):
    """
    Multiplies A by its Hermitian conjugate on last two indices

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    ndarray
        Same shape as input array.

    """
    return np.matmul(A, hconj(A))

def mult_trans(A):
    """
    Multiplies A by its transpose on last two indices

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    ndarray
        Same shape as input array.

    """
    return np.matmul(A, trans(A))

def mat2vec(A):
    sz = A.size 
    return A.reshape(sz)

def vec2mat(Av):
    sz = Av.size 
    return Av.reshape((sz//9, 3, 3))

def eig_vals(A):
    """
    Extracts eigenvalues of array A

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    ndarray
        shape (m,n,...,3).

    """
    A_dim = A.ndim
    A_shape = A.shape
    A_size = A.size
    A_flatter = A.reshape(A_size//9, 3, 3)
    e_vals = np.zeros((A_size//9, 3), dtype=complex)
    for n, A_el in enumerate(A_flatter):
        e_vals[n,:] = sl.eigvals(A_el)
    return e_vals.reshape(A_shape[0:A_dim-2] + (3,))


def eig_orbital(A):
    """
    Extracts eigenvalues of array A^dagger A, relevant for angular momentum vectors

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    ndarray
        shape (m,n,...,3).

    """
    H = np.matmul(hconj(A), A)
    H_dim = H.ndim
    H_shape = H.shape
    H_size = H.size
    H_flatter = H.reshape(H_size//9, 3, 3)
    e_vals = np.zeros((H_size//9, 3), dtype=float)
    for n, H_el in enumerate(H_flatter):
        e_vals[n,:] = sl.eigvals(H_el)
    e_vals.sort()
    return e_vals.reshape(H_shape[0:H_dim-2] + (3,))


def eig_spin(A):
    """
    Extracts eigenvalues of array A A^dagger, relevant for spin vectors

    Parameters
    ----------
    A : ndarray dtype complex, shape (m,n,...,3,3)
        order parameter array.

    Returns
    -------
    ndarray
        shape (m,n,...,3).

    """
    H = np.matmul(A, hconj(A))
    H_dim = H.ndim
    H_shape = H.shape
    H_size = H.size
    H_flatter = H.reshape(H_size//9, 3, 3)
    e_vals = np.zeros((H_size//9, 3), dtype=float)
    for n, H_el in enumerate(H_flatter):
        e_vals[n,:] = sl.eigvals(H_el)
    e_vals.sort()
    return e_vals.reshape(H_shape[0:H_dim-2] + (3,))

def norm(D, n=1):
    """
    n-Norm of complex square array D, defined as $tr(D D^\dagger)^{n/2}$

    Parameters
    ----------
    D : Complex array shape 
        Order parameter array.

    Returns
    -------
    float
        n-Norm of array.

    """
    # D_H = np.conj(np.swapaxes(D, dim-2, dim-1))
    D_H = hconj(D)
    norm2 = tr(np.matmul(D, D_H)).real
    # if n == 1:
    #     return norm2**(n/2)
    # else:
    return norm2**(n/2)


def inner(X, Y):
    """
    Inner product of complex square arrays X, Y

    Parameters
    ----------
    X,Y : Complex array shape 
        Order parameter array.

    Returns
    -------
    float, dtype complex
        Inner product.

    """
    dim = Y.ndim
    # X_H = np.conj(np.swapaxes(Y, dim-2, dim-1))
    Y_H = np.conj(np.swapaxes(Y, dim-2, dim-1))
    return (tr(np.matmul(X, Y_H))).real


def distance(X, Y):
    """
    Distance between (3,3) complex arrays defined by norm.    

    Parameters
    ----------
    X : Complex array shape (3,3)
        Order parameter array.
    Y : Complex array shape (3,3)
        Order parameter array.

    Returns
    -------
    float
        Distance between X and Y.

    """
    D = X - Y
    return norm(D)


def project(X, Y):
    """
    Projects X onto Y, X, Y (3,3) complex arrays.
    X - (X,Y) Y/(Y,Y)
       

    Parameters
    ----------
    X : Complex array shape (3,3)
        Order parameter array.
    Y : Complex array shape (3,3)
        Order parameter array.

    Returns
    -------
    float
        X projected onto Y.

    """

    return X -  np.multiply.outer(inner(X,Y)/inner(Y,Y), Y )


def lengths(X, Y):
    """
    Length of X along Y and orthoginal to Y
    Projects X onto Y, X, Y (3,3) complex arrays.
    X - (X,Y) Y/(Y,Y)
       

    Parameters
    ----------
    X : Complex array shape (n,m,..3,3)
        Order parameter array.
    Y : Complex array shape (3,3)
        Order parameter array.

    Returns
    -------
    float
        X projected onto Y.

    """
    print((inner(X,Y)/inner(Y,Y)).shape)
    print(Y.shape)
    X_Y = X - np.multiply.outer(inner(X,Y)/inner(Y,Y), Y ) 

    return np.array([norm(X_Y), norm(X - X_Y)]).T

def levi_civita3_matrix2vector(X):
    X_dim = X.ndim
    X_shape = X.shape
    X_size = X.size
    X_flatter = X.reshape(X_size//9, 3, 3)
    eps_X = np.array([X_flatter[:,1,2] - X_flatter[:,2,1], 
                      X_flatter[:,2,0] - X_flatter[:,0,2], 
                      X_flatter[:,0,1] - X_flatter[:,1,0]])
    return eps_X.reshape(X_shape[0:X_dim-2] + (3,)) 

def levi_civita3_vector2matrix(v):
    v_dim = v.ndim
    v_shape = v.shape
    v_size = v.size
    v_flatter = v.reshape(v_size//3, 3)
    eps_X = np.zeros((v_size//3, 3, 3))
    eps_X[:,0,1] = v_flatter[:,2]
    eps_X[:,1,2] = v_flatter[:,0]
    eps_X[:,0,2] = -v_flatter[:,1]
    eps_X[:,1,0] = -v_flatter[:,2]
    eps_X[:,2,1] = -v_flatter[:,0]
    eps_X[:,2,0] = +v_flatter[:,1]

    return eps_X.reshape(v_shape[0:v_dim-1] + (3,3)) 