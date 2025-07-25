# -*- coding: utf-8 -*-

import numpy as np
import he3_props as h3p
import he3_bases as h3b
import he3_constants as h3c
import he3_data as h3d
import he3_matrix as h3m
import he3_free_energy as h3f


# Theory functions
def s_scale(p):

    """ Natural entropy density scale, units joule per kelvin per nm3 .
    """
    return h3p.f_scale(p) / (h3p.T_mK(1, p) * 1e-3)


def entropy_density_norm(t, p, phase):
    """
    Normalised entropy density, in units of $(1/3)N(0)k_B^2T_c$.

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    p : TYPE
        DESCRIPTION.
    phase : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    A = h3p.delta_phase_norm(t, p, phase) * h3b.D_dict[phase]
    dim = A.ndim
    A_T = np.swapaxes(A, dim-2, dim-1)
    # A_C = np.conj(A)
    A_H = np.conj(A_T)

    AA_H = np.matmul(A , A_H)
    AA_T = np.matmul(A , A_T)
        
    delta_beta_arr = h3c.beta_const * h3p.delta_beta_norm_asarray(p)
    
    dUn0_dT = h3m.tr( AA_H )
    
    dUn1_dT = delta_beta_arr[0] *  h3m.tr( AA_T) * np.conj(h3m.tr( AA_T))
    dUn2_dT = delta_beta_arr[1] *  h3m.tr( AA_H )**2
    dUn3_dT = delta_beta_arr[2] *  h3m.tr( np.matmul( AA_T , np.conj(AA_T))  )
    dUn4_dT = delta_beta_arr[3] *  h3m.tr( np.matmul( AA_H , AA_H )  )
    dUn5_dT = delta_beta_arr[4] *  h3m.tr( np.matmul( AA_H , np.conj(AA_H) )  )
    
    # This is entropy density relative to the normal phase, so need to add it 
    # Normal phase value is $\pi^2 t$ in these units.
    
    s_norm = np.pi**2 * t + dUn0_dT + dUn1_dT + dUn2_dT + dUn3_dT + dUn4_dT + dUn5_dT
    
    return s_norm.real

def specific_heat_N_norm(t, p):
    """
    Specific heat in the normal phase, in of $(1/3)N(0)k_B^2T_c$.

    Parameters
    ----------
    t : TYPE
        Reduced temperature, $T/T_c$.
    p : TYPE
        Pressure in bar.

    Returns
    -------
    Specific heat, same type as t.

    """
    
    return np.pi**2 * t

def enthalpy_density_norm(t, p, phase):
    
    return t * entropy_density_norm(t, p, phase)

def energy_density_norm(t, p, phase):
    
    return enthalpy_density_norm(t, p, phase) - h3p.f_phase_norm(t, p, phase)
    
def theta_norm(t, p, phase):
    
    return 0.25*enthalpy_density_norm(t, p, phase) - h3p.f_phase_norm(t, p, phase)

def alpha_pt_AB_norm(t, p):

    delta_theta = theta_norm(t, p, 'A') - theta_norm(t, p, 'B')
    
    thermal_energy_density = 0.75 * enthalpy_density_norm(t, p, 'A')
    
    return delta_theta/thermal_energy_density

