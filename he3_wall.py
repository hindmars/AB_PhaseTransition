#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:20:44 2021

@author: hindmars
"""

"""
Simple routines for finding the critical bubble in He3 potential. 
Uses newton-krylov solver.

Default potential is scaled so that: 
    stable (broken) phase has order parameter phi = 1, 
    stable (broken) phase has potential = 0
    potential grows phi**4/4

Then there is only one free parameter, the value of phi at the maximum, phi_m, 
and

0 < phi_m < 0.5

As phi_m -> 0.5 we approach the thin wall limit and the bubble grows 
arbitraily large.

As phi_m -> the barrier disappears and the critical bubble becomes a small bump of 
vanishing low energy.

To get started, try:
    import critical_bubble_minimal as cbm
    phi, pot, gr = cbm.krylov_bubble(0.4,display=True)

First argument is phi_m.

The computation is carried out on a default grid with 200 points and max radius 20
You can change these values with e.g.

    phi, pot, gr = cbm.krylov_bubble(0.45, gr_pars=(400,40))

(this is more like a thin wall bubble)

It also works in dimensions 1, 2, 3 and 4, e.g.

    phi, pot, gr = cbm.krylov_bubble(0.45, gr_pars=(400,40), dim=4)

If you want to know the energy

    cbm.energy(phi, pot, gr)

Enjoy!

Mark
3.6.21

MH 13.7.21: One can get initialise the potential with a more conventional form of the potential, 
resembling the Ginzburg-Landau free energy, by creating a quartic_paramset object, 
and invoking the method  quartic_params_init. See accompanying file critical_bubble_demo

"""


import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import newton_krylov
from scipy.integrate import solve_ivp
import scipy.optimize.nonlin

import he3_tools as  h

Kxx = np.diag([3,1,1]) 

class material_paramset:
    """
    Class for the material parameters of the static GL functional, including temperature and pressure dependence.

    A material_paramset object can also be directly passed to quartic_potential, along with a temperature.
    """
    def __init__(self, *args):
        al, bn, t, p = h.args_parse(*args)
        self.alpha = al
        self.beta_arr = bn
        self.t = t
        self.p = p

    def beta_A_norm(self):
        return np.sum( self.beta_arr * h.R_dict["A"] )
        
    def beta_B_norm(self):
        return np.sum( self.beta_arr * h.R_dict["B"] )

    def beta_phase_norm(self, phase):
        return np.sum( self.beta_arr * h.R_dict[phase] )

    def delta_A_norm(self):
        return np.sqrt(- self.alpha/(2 * self.beta_A_norm()))
        
    def delta_B_norm(self):
        return np.sqrt(- self.alpha/(2 * self.beta_B_norm()))

    def delta_phase_norm(self, phase):
        return np.sqrt(- self.alpha/(2 * self.beta_phase_norm(phase)))

    def delta_wc(self):
        return np.sqrt(- self.alpha/(2 * h.beta_const))

    def f_A_norm(self):
        return -0.25* self.alpha**2 /( self.beta_A_norm())
        
    def f_B_norm(self):
        return -0.25* self.alpha**2 /( self.beta_B_norm())

    def f_phase_norm(self, phase):
        return -0.25* self.alpha**2 /( self.beta_phase_norm(phase))

    def xi(self):
        """Ginzburg Landau correlation length.
        """
        if self.p is not None:
            xiGL = h.xi(self.t, self.p)
        else:
            xiGL = h.xiGL_const/np.sqrt(-self.alpha)
        return xiGL


class quartic_potential:
    """ 
    Quartic potential class. 

    """
    
    def __init__(self, *args):
        # self.t = t
        # self.p = p
        self.args = args
        # al, bn = h.args_parse(args)
        # self.alpha = al
        # self.beta_arr = bn
        self.mat_pars = material_paramset(*args)

    def v(self, A):
        return h.U(A, self.mat_pars.alpha, self.mat_pars.beta_arr)
    
    def v_prime(self, A):
        return h.dU_dA(A, self.mat_pars.alpha, self.mat_pars.beta_arr)
        
    # def plot_line_section(self, X, Y, scale=None, n=500):
    #     v, A_XD, U_XD = h.line_section(X, Y, self.args, scale, n)
    #     fig, ax = plt.subplots()
    #     ax.plot(v, U_XD)
    #     ax.set_xlabel(r'$\phi/\Delta_{\rm wc}$')
    #     ax.set_ylabel(r'$U(\phi)$')
    #     ax.grid()
    #     return ax

    
class grid_1d:
    """ 1D grid class. Contains n x values betwen 0 and R, at half dx intervals. 
    Also values x_minus (shifted grid).  dim is number of space dimensions 
    in differential equation, not grid dimension. 
    Thermal activation in d=3: dim = 3
    Quantum tunnelling in d=4: dim = 4
    """
    
    def __init__(self, n, R, dim=3):
        self.n = n
        self.R = R
        self.dx = R/n # unidorm grid for now.
        self.dx2 = self.dx**2
        self.dim = dim
    
        x_min = 0.5*self.dx
        x_max = R - x_min
    
        self.x = np.linspace(x_min, x_max, n) # uniform grid for now.
    
        self.x_minus = np.roll(self.x,1)
        self.x_minus[0] = self.x[0]


class thin_wall_bubble:
    
    def __init__(self, pot, dim=3):
        self.surface_energy = np.sqrt(h.beta_const) * pot.mat_pars.delta_wc()**3/6
        self.energy_diff = pot.mat_pars.f_A_norm() - pot.mat_pars.f_B_norm()
        if self.energy_diff > 0:
            if dim == 4:
                self.r_bub = 3*self.surface_energy/self.energy_diff
                self.e_bub = (27*np.pi**2/2) * self.surface_energy**4/self.energy_diff**3
            if dim == 3:
                self.r_bub = 2*self.surface_energy/self.energy_diff
                self.e_bub = (16*np.pi/3) * self.surface_energy**3/self.energy_diff**2
            elif dim == 2:
                self.r_bub = self.surface_energy/self.energy_diff
                self.e_bub = np.pi * self.surface_energy**2/self.energy_diff
            elif dim == 1:
                self.r_bub = np.inf
                self.e_bub = np.inf
        else:
            self.r_bub = np.inf
            self.e_bub = np.inf


def field_eqn(phi, pot, gr):
    """
    Scalar field equation in 1 radial dimension, with trivial gradient term.
    """
    
    phi_plus = np.roll(phi,-1, axis=0)
    phi_minus = np.roll(phi,+1, axis=0)

    # Neumann boundary conditions - derivative vanishes at r=0, max(r)    
    phi_plus[-1] = phi[-2]
    phi_minus[0] = phi[1]
    
    # Need to sort out gradient term

    second = (phi_plus + phi_minus - 2*phi)/gr.dx2
    first = np.matmul((phi_plus  - phi_minus) , np.multiply.outer(1/(gr.x + gr.x_minus), h.id3))
    first *= (gr.dim - 1)/gr.dx
    
    return np.matmul(Kxx, second + first) - pot.v_prime(phi)


def initial_condition(pot, gr, D=None):
    """
    Initial guess: a kink-antikink pair symmetric around r=0, at +/- r_bub
    """
    
    bub = thin_wall_bubble(pot, dim=gr.dim)
    
    m2 = 1/pot.mat_pars.xi()
    k = np.sqrt(m2/4)

    rb = bub.r_bub
    if rb == np.inf:
        rb = gr.R/2

    phi_init = 0.25*(1 - np.tanh(k*(gr.x - rb)))*(1 + np.tanh(k*(gr.x + rb )))  

    if D is None:
        D = pot.mat_pars.delta_A_norm()*h.D_dict["A"] - pot.mat_pars.delta_B_norm()*h.D_dict["B"]
    # D = D/h.norm(D)

    if gr.dim == 1:
        A_init = h.D_dict["B"] * pot.mat_pars.delta_B_norm() \
            + np.multiply.outer(phi_init , D ) 
    else:
        A_init = h.D_dict["A"] * pot.mat_pars.delta_A_norm() \
            - np.multiply.outer(phi_init , D ) 
        
        
    return A_init


def krylov_bubble(*args, gr_pars=(200,20), dim=3, display=False):
    """
    Apply Krylov solver to find unstable bubble solution. 
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    

    # t = param[0]
    # p = param[1]        
    
    pot = quartic_potential(*args)
    
    gr = grid_1d(*gr_pars, dim=dim)
    # bub = thin_wall_bubble(pot, dim=dim)
    
    phi_init = initial_condition(pot, gr)
    
    def field_eqn_fix(phi):
        return field_eqn(phi, pot, gr)

    try:
        phi = newton_krylov(field_eqn_fix, phi_init, verbose=True, maxiter=200)
    except scipy.optimize.nonlin.NoConvergence as e:
        phi = e.args[0]
        print('No Convergence')

    # if display:
        
    #     plt.figure()
    #     plt.plot(gr.x, phi_init/h.delta_wc, label='Thin wall') 
    #     plt.plot(gr.x, phi/h.delta_wc, label='Final' )
    #     plt.xlabel(r'$r$')
    #     plt.ylabel(r'$\phi/\phi_{\rm b}$')
    #     plt.xlim(0,gr.R)
    #     plt.grid()
    #     plt.legend()

    return phi, pot, gr 

def relax(t_eval, *args, gr_pars=(200,20), dim=1):
    """
    Apply relaxation method to find domain wall solution. 
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    
    
    pot = quartic_potential(*args)
    
    gr = grid_1d(*gr_pars, dim=dim)
    # bub = thin_wall_bubble(pot, dim=dim)
    
    phi_init = initial_condition(pot, gr)
    
    def field_eqn_vec(tau, phi_vec):
        force_mat = field_eqn(h.vec2mat(phi_vec), pot, gr)
        return h.mat2vec(force_mat)

    tspan = [min(t_eval), max(t_eval)]
    sol = solve_ivp(field_eqn_vec, tspan, h.mat2vec(phi_init), t_eval=t_eval)

    Asol = sol.y
    print(Asol.shape)
    return Asol.reshape(gr.n, 3, 3, len(t_eval)), pot, gr 

def energy_density(phi, pot, gr):
    phi_plus = np.roll(phi,-1, axis=0)
    phi_plus[-1] = phi[-2]
    
    first = (phi_plus  - phi)/gr.dx
    # Need to sort out gradient term
    eden_grad = h.tr(np.dot(first, Kxx) @ h.hconj(first) ).real
    eden_pot = pot.v(phi)-pot.v(phi[-1])
    
    return eden_grad + eden_pot, eden_grad, eden_pot

def energy(phi, pot, gr):

    if gr.dim == 3:
        vol_factor = 4*np.pi/3
    elif gr.dim == 4:
        vol_factor = 2*np.pi**2
    elif gr.dim == 1:
        vol_factor=1    
    elif gr.dim == 2:
        vol_factor=2*np.pi    
    else:
        raise ValueError("energy: only with grid.dim = 1,2,3,4 at the moment")

    eden_tot, eden_grad, eden_pot = energy_density(phi, pot, gr)
    
    e_grad = vol_factor * np.trapz(eden_grad, gr.x**(gr.dim) )
    e_pot = vol_factor * np.trapz(eden_pot, gr.x**(gr.dim))
    
    return e_grad + e_pot, e_grad, e_pot


def action(phi, pot, gr):
    
    e_tot, e_grad, e_pot = energy(phi, pot, gr)
    
    if gr.dim == 3:
        if pot.T is not None:
            scale = pot.T
        else:
            scale = pot.phi_b
            print("No temperature defined - dividing bubble energy by phi_b instead of T$")
    elif gr.dim == 4:
        scale = 1.0
    else:
        raise ValueError("action: only with grid.dim = 3,4 at the moment")

    return e_tot/scale, e_grad/scale, e_pot/scale


