#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:20:44 2021

@author: hindmars; modified by timohyva@github
"""

"""
Simple routines for finding the critical bubble in He3 potential. 
Uses newton-krylov solver.

Default potential is scaled so that: 
    stable (broken) phase has order parameter phi = 1 i.e., true vacuum (superfluid B-Phase)
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
import matplotlib.pyplot as plt
from scipy.optimize import newton_krylov
from scipy.integrate import solve_ivp
import scipy.optimize.nonlin
import scipy.linalg as sl

import he3_tools as  h

# Stiffness parameter for configurations dependinng only on x
Kxx = np.diag([3,1,1]) 

# Various standard boundary conditions
bc_dir = [np.array([1,1,1]), np.array([0,0,0])]
bc_neu = [np.array([0,0,0]), np.array([1,1,1])]
bc_rob = [np.array([1,1,1]), np.array([1,1,1])]
bc_min_pb = [np.array([1,0,0]), np.array([0,1,1])]
bc_med_pb = [np.array([1,1,1]), -np.array([0,1,1])]
bc_max_pb = bc_dir
# Diffuse bcs from Wiman and Sauls 2016, reporting Ambegaokar, de Gennes, Rainer 1975
b_adgr = 0.54
bc_diff = [np.array([1, h.xiGL_const, h.xiGL_const]), -np.array([0, b_adgr, b_adgr])]
bcleft = 0
bcright = -1

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

    def xi0(self):
        """Ginzburg Landau correlation length.
        """
        if self.p is not None:
            xiGL = h.xi(0, self.p)
        else:
            xiGL = h.xiGL_const
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
    
    def __init__(self, n, R, dim=3, bcs=[bc_neu, bc_neu]):
        self.n = n
        self.R = R
        self.dx = R/n # unidorm grid for now.
        self.dx2 = self.dx**2
        self.dim = dim
        self.bcs = bcs
    
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
    if gr.dim > 1:
        first = np.matmul((phi_plus  - phi_minus) , 
                          np.multiply.outer(1/(gr.x + gr.x_minus), h.id3))
        first *= (gr.dim - 1)/gr.dx
        derivs = second + first
    else:
        derivs = second
        
    return np.matmul(Kxx, derivs) - pot.v_prime(phi)

def find_boundary_condition(bcs, boundary):
    a1 = bcs[boundary][0]
    a2 = bcs[boundary][1]
    dirichlet = (a2 == 0)
    neumann = (a1 == 0)
    robin = np.logical_not(np.logical_or(dirichlet, neumann))
    return dirichlet, neumann, robin    

# def boundary_condition(A, bcs, boundary, grid):
def boundary_condition(A, grid, boundary):
    
    bcs = grid.bcs
    dirichlet, neumann, robin = find_boundary_condition(bcs, boundary)

    if boundary == bcleft:
        neighbour = 1
    elif boundary == bcright:
        neighbour = -2
    else:
        raise ValueError("boundary must be 0 (left) or -1 (right)")

    if np.any(dirichlet):
        A[boundary, :, dirichlet] = 0
    
    # Following should be incorporated into shifted phi, as they involve derivatives?
    # if np.any(neumann):
    #     A[boundary, :, neumann] = A[neighbour, :, neumann] 
    # if np.any(robin):
    #     a1 = bcs[boundary][0]
    #     a2 = bcs[boundary][1]
    #     A[boundary, :, robin] = A[neighbour, :, robin] / (1 - grid.dx * (a1[robin]/a2[robin]))
    
    return A

def phi_shift(phi, direction, gr):
    
    bcs = gr.bcs
        
    phi_shifted = np.roll(phi,-direction, axis=0)
    
    if direction == +1:
        boundary = bcright
    elif direction == -1:
        boundary = bcleft
    else:
        raise ValueError('phi_shift direction must be +1 or -1')
    
    a1 = bcs[boundary][0]
    a2 = bcs[boundary][1]
    dirichlet = (a2 == 0)
    neumann = (a1 == 0)
    robin = np.logical_not(np.logical_or(dirichlet, neumann))

    if boundary == bcleft:
        neighbour = 1
    elif boundary == bcright:
        neighbour = -2
    else:
        raise ValueError("boundary must be 0 (left) or -1 (right)")

    if np.any(dirichlet):
        # phi_shift[boundary, :, dirichlet] = phi[boundary, :, dirichlet] - phi[neighbour, :, dirichlet]
        # Assumes Dirichlet BC is phi = 0
        phi_shifted[boundary, :, dirichlet] =  - phi[neighbour, :, dirichlet]
    if np.any(neumann):
        phi_shifted[boundary, :, neumann] = phi[neighbour, :, neumann]
    if np.any(robin):
        a12 = a1[robin]/a2[robin]
        factor_shape_21 = ((1 + a12*gr.dx)/(1 - a12*gr.dx))[:,None]
        phi_shifted[boundary, :, robin] = phi[neighbour, :, robin]*factor_shape_21

    return phi_shifted


def field_eqn_with_bcs(phi, pot, gr):
    """
    GL field equation in 1 radial dimension
    """

    # phi_plus = np.roll(phi,-1, axis=0)
    # phi_minus = np.roll(phi,+1, axis=0)
    # phi_plus[-1] = phi[-2]
    # phi_minus[0] = phi[1]
    phi_plus = phi_shift(phi, +1, gr)
    phi_minus = phi_shift(phi, -1, gr)
    
    phi = boundary_condition(phi, gr, bcleft)
    phi = boundary_condition(phi, gr, bcright)
    # phi_plus = phi_shift(phi, +1, bcs)
    # phi_minus = phi_shift(phi, -1, bcs)
    # else:
    #     phi_plus = np.roll(phi,-1, axis=0)
    #     phi_minus = np.roll(phi,+1, axis=0)
    #     phi_plus[-1] = phi[-2]
    #     phi_minus[0] = phi[1]
        
    second = (phi_plus + phi_minus - 2*phi)/gr.dx2 

    if gr.dim > 1:
        first = np.matmul((phi_plus  - phi_minus) , 
                          np.multiply.outer(1/(gr.x + gr.x_minus), h.id3))
        first *= (gr.dim - 1)/gr.dx
        derivs = second + first
    else:
        derivs = second
        
    return np.matmul(Kxx, derivs) - pot.v_prime(phi)



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

def expm_vec(x, T):
    """ Calculates exp(x * T), where T is a matrix and x an array.
    Returns array with shape x.shape + T.shape.
    """
    foo = np.zeros(x.shape + T.shape, dtype=complex)
    for n,x in enumerate(x):
        foo[n,:,:] = sl.expm(x * T )
    return foo


def apply_texture(A, gr, T_list, fun_list=[expm_vec]*3):
    """ Applies textture by applying O(3) rotations and phase angle.
    T_list[2] should be 1j*h.id3. 
    signature of 3 functions in funlist should be (gr.x, T)
    """
    R_spin  = fun_list[0](2*np.pi*gr.x/gr.R, T_list[0])
    R_angm  = fun_list[1](2*np.pi*gr.x/gr.R, T_list[1])
    R_phase = fun_list[2](2*np.pi*gr.x/gr.R, T_list[2])
    
    return np.matmul(np.matmul(np.matmul(R_spin, A), R_angm), R_phase)


def initial_condition_confine(pot, gr, wall_phase="Ay", bulk_phase="B", 
                              # bcs=[bc_dir, bc_dir],
                              T_list=None,
                              fun_list=[expm_vec]*3):
    """
    Initial guess for order parameter. Interpolates between wall and bulk 
    phases over GL coherence length.  If RHS b.c. is Neumann, it is left in bulk phase.
    Applies texture if T_list is none
    """
    bcs = gr.bcs
    # bub = thin_wall_bubble(pot, dim=gr.dim)
    
    m = pot.mat_pars.xi0()/pot.mat_pars.xi()
    k = m/2

    w = gr.R
    phi_init_left = np.tanh(k*gr.x)
    phi_init_right = np.tanh(k*(w - gr.x))

    phi_init = np.ones_like(gr.x)
    proj = np.multiply.outer(np.ones_like(phi_init), h.id3)
    # Test to see if we have dirichlet at boundaries, in which case we apply tanh profile there
    
    dch_l, _, _ = find_boundary_condition(bcs, bcleft)
    dch_r, _, _ = find_boundary_condition(bcs, bcright)
    
    if np.any(dch_l):
        print("Initial condition: Applyng Dir to left boundary")
        phi_init *= phi_init_left
        dch_l_vec = dch_l.astype(float)
        not_dch_l_vec = np.logical_not(dch_l).astype(float)
        proj_l = np.multiply.outer(phi_init_left, np.diag(dch_l_vec) ) + \
            np.multiply.outer(np.ones_like(phi_init_left), np.diag(not_dch_l_vec) )  
        proj = np.matmul(proj_l, proj)
    if np.any(dch_r):
        print("Initial condition: Applyng Dir to right boundary")
        phi_init *= phi_init_right
        dch_r_vec = dch_r.astype(float)
        not_dch_r_vec = np.logical_not(dch_r).astype(float)
        proj_r = np.multiply.outer(phi_init_right, np.diag(dch_r_vec) ) + \
            np.multiply.outer(np.ones_like(phi_init_right), np.diag(not_dch_r_vec) )  
        proj = np.matmul(proj_r, proj)

    D = pot.mat_pars.delta_phase_norm(bulk_phase)*h.D_dict[bulk_phase] - \
            pot.mat_pars.delta_phase_norm(wall_phase)*h.D_dict[wall_phase]
    A_init = h.D_dict[wall_phase] * pot.mat_pars.delta_phase_norm(wall_phase) \
            + np.multiply.outer(phi_init , D ) 

    # Minimal pair-breaking inital condition (Dirichlet for normal component)
    # proj = np.multiply.outer(phi_init, np.diag([1,0,0]) ) + \
    #     np.multiply.outer(np.ones_like(phi_init), np.diag([0,1,1]) ) 

    A_init = np.matmul(proj, A_init)
    
    if T_list is not None:
        A_init = apply_texture(A_init, gr, T_list, fun_list)
    
    return A_init

def krylov_bubble(*args, gr_pars=(200,20), dim=3, 
                  display=False, 
                  verbose=True, maxiter=200, **kwargs):
    """
    Apply Krylov solver to find unstable bubble or wall (dim=1) solution. 
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    

    pot = quartic_potential(*args)
    
    gr = grid_1d(*gr_pars, dim=dim, bcs=[bc_neu, bc_neu])
    
    phi_init = initial_condition(pot, gr)

    
    def field_eqn_fix(phi):
        return field_eqn(phi, pot, gr)

    try:

        phi = newton_krylov(field_eqn_fix, phi_init, verbose=verbose, maxiter=maxiter, **kwargs)
    except scipy.optimize.nonlin.NoConvergence as e:
        phi = e.args[0]
        print('No Convergence')

    return phi, pot, gr 


def krylov_confine(*args, gr_pars=(200,20), dim=1, 
                   wall_phase="Ay", bulk_phase="B", 
                   bcs = [bc_min_pb, bc_min_pb], 
                   T_list = None, **kwargs):
    """
    Apply Krylov solver to find solution in confined boundary conditions.
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    
    pot = quartic_potential(*args)
    gr = grid_1d(*gr_pars, dim=dim, bcs=bcs)
    
    phi_init = initial_condition_confine(pot, gr, 
                                         wall_phase=wall_phase, bulk_phase=bulk_phase,
                                         # bcs=bcs, 
                                         T_list=T_list)
    
    def field_eqn_fix(phi):
        return field_eqn_with_bcs(phi, pot, gr)

    try:
        phi = newton_krylov(field_eqn_fix, phi_init, **kwargs)
    except scipy.optimize.nonlin.NoConvergence as e:
        phi = e.args[0]
        print('No Convergence')

    return phi, pot, gr 


def relax(t_eval, *args, gr_pars=(200,20), dim=1):

    """
    Apply relaxation method to find domain wall solution. 
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    
    
    pot = quartic_potential(*args)
    
    gr = grid_1d(*gr_pars, dim=dim)
   
    phi_init = initial_condition(pot, gr)

    
    def field_eqn_vec(tau, phi_vec):
        force_mat = field_eqn(h.vec2mat(phi_vec), pot, gr)
        return h.mat2vec(force_mat)


    tspan = [min(t_eval), max(t_eval)]
    sol = solve_ivp(field_eqn_vec, tspan, h.mat2vec(phi_init), t_eval=t_eval)

    Asol = sol.y
    print(Asol.shape)
    return Asol.reshape(gr.n, 3, 3, len(t_eval)), pot, gr 



def relax_from_ic(t_eval, phi_init, pot, gr):
    """
    Apply relaxation method to find domain wall solution. 
    
    Returns: order parameter phi with extra t dimension, the potential object, and the grid object.
    """

    def field_eqn_vec(tau, phi_vec):
        force_mat = field_eqn_with_bcs(h.vec2mat(phi_vec), pot, gr)
        return h.mat2vec(force_mat)

    tspan = [min(t_eval), max(t_eval)]
    sol = solve_ivp(field_eqn_vec, tspan, h.mat2vec(phi_init), t_eval=t_eval)

    Asol = sol.y
    print(Asol.shape)
    Asol = Asol.reshape(gr.n, 3, 3, len(t_eval))
    print(Asol.shape)
    Apg_list = []
    for n,t in enumerate(t_eval):
        Apg_list.append((Asol[:,:,:,n], pot, gr))
    return Apg_list


def energy_density(phi, pot, gr):
    phi_plus = np.roll(phi,-1, axis=0)
    phi_plus[-1] = phi[-2]
    
    first = (phi_plus  - phi)/gr.dx
    # Need to sort out gradient term
    eden_grad = h.tr(np.dot(first, Kxx) @ h.hconj(first) ).real
    eden_pot = pot.v(phi)
    
    return eden_grad + eden_pot, eden_grad, eden_pot

def surface_energy(A, pot, gr):
    """Calculates surface energy of configuration, relative to minimum energy density.
    """
    try:
        bcs = gr.bcs
    except:
        bcs = [bc_neu, bc_neu]    
    
    eden, eden_grad, eden_pot = energy_density(A, pot, gr)
    en = np.trapz(eden - np.min(eden), gr.x)
    
    if bcs is None:
        n_surface = 1
    else:
        dch, neu, rob = find_boundary_condition(bcs, bcright)
    
        if np.all(neu):
            n_surface = 1
        else:
            n_surface = 2

    # print("Surface energy called with n_surface =", n_surface)

    return en/n_surface 

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
    
    return np.array([e_grad + e_pot, e_grad, e_pot])


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


def plot_eigs(A, pot, gr, angm="orbital"):
    
    t = pot.mat_pars.t
    p = pot.mat_pars.p
    eden, eden_grad, eden_pot = energy_density(A, pot, gr)
    # sigma_tot = np.trapz(eden, gr.x)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    sigma_bw = surface_energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))

    fig, ax = plt.subplots(2,1, sharex='col')
    # gridspec_kw={'hspace': 0, 'wspace': 0}
    x =  gr.x * h.xi(0,p)/h.xi(t,p)
    xmin = min(x)
    xmax = max(x)
    
    w_mu = gr.R * h.xi(0,p)/1000
    
    ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm()) + 1)
    
    ax[0].set_ylabel(r'$e/|f_B| + 1$')
    ax[0].grid()
    ax[0].set_xlim(xmin, xmax)
    tstring1 = r'p={:.1f} bar, $T={:.2f}$ mK, $w={:.2f} \;\mu$m'.format(p, t*h.Tc_mK(p), w_mu)
    tstring2 = r'$\sigma_{{bw}}/\xi_{{\rm GL}}(T)|f_B(T)| = {:.3f}$, '.format(sigma_bw)
    
    #Now plot eigenvalues
    if angm=="orbital":
        eigs = h.eig_orbital(A)
        tstring3 = r'$A^\dagger A v_a = \lambda_a v_a}$, orbital'
    elif angm=="spin":
        eigs = h.eig_spin(A)
        tstring3 = r'$A A^\dagger v_a = \lambda_a v_a}$, spin'
    else:
        raise ValueError("angm must be orbital or spin")
        
    ax[0].set_title(tstring1 + '\n' + tstring2 + tstring3 )
    norm =  np.sqrt(3)/h.delta_B_norm(t, p) 
    
    ax[1].plot(x, eigs[:,0]/norm**2, label=r'$\lambda_1$')
    ax[1].plot(x, eigs[:,1]/norm**2, label=r'$\lambda_2$')
    ax[1].plot(x, eigs[:,2]/norm**2, label=r'$\lambda_3$')
    
    ax[1].set_xlim(xmin, xmax)
    ax[1].legend(loc='best')
    ax[1].grid()
    ax[1].set_ylabel(r'$3\lambda_a/\Delta_B^2(T,p)$')
    ax[1].set_xlabel(r'$x/\xi_{\rm GL}(T)$')

    return ax


def thuneberg_formula(t,p):
    """Formula for AB interface free energy from Thuneberg PRB 1991.
    """
    beta1 = h.beta_norm(t,p,1)
    beta2 = h.beta_norm(t,p,2)
    beta3 = h.beta_norm(t,p,3)
    beta45 = h.beta_norm(t,p,4) + h.beta_norm(t,p,5)
    beta345 = beta3 + beta45

    beta2_0 = beta1 + beta2 + beta345/2
    
    beta3_0 = beta2_0
    beta1_0 = -beta2_0/2
    beta34_0 = 2*beta2_0
    
    a = 2*beta1 + beta2 - beta45
    c = -(2*beta1 + beta345)
    
    kappa = (beta3_0*(beta1_0 + 3*beta2_0)/(4*beta2_0*beta34_0))**0.5
    I2 = 1.89 - 1.98*kappa**0.5 - 0.31*kappa
    
    if c > 0:
        I1 = (a+c)**0.5 + (a/c**0.5)*np.log(((a+c)**0.5 + c**0.5)/a**0.5)
    elif c < 0:
        I1 = (a+c)**0.5 + (a/(-c)**0.5)*np.asin(((-c)**0.5)/a**0.5)
    else:
        I1 = np.nan
    
    
    square_bracket_20 = I1/(2*beta2_0)**0.5 + \
        0.5*I2*(4*a**2/(beta2_0*beta3_0*(beta1_0 + beta2_0)))**0.25
        
    return square_bracket_20*h.xi(t,p) * h.alpha_norm(t)**2/(4*beta2_0)

def get_wall(t, p, w, N=500, **kwargs):
    
    L = w/h.xi(0,p)
    # bcs = [bc_neu, bc_neu]
    
    A, pot, gr = krylov_bubble(t, p, gr_pars=(N,L), dim=1, verbose=False, **kwargs)

    sigma_bw = surface_energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    sigma_tot = energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    # sigma_bw = sigma_bw/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))
    
    return (A, pot, gr), sigma_bw, sigma_tot

def plot_wall(A, pot, gr, 
              real_comp_list=[(0,0), (1,1), (2,2)], 
              imag_comp_list = [(0,1), (1,0)],
              legend_loc='center right'):
    comp_key = ('x', 'y', 'z')
    t = pot.mat_pars.t
    p = pot.mat_pars.p
    eden, eden_grad, eden_pot = energy_density(A, pot, gr)
    sigma = surface_energy(A, pot, gr)*h.xi(0,p)/(abs(pot.mat_pars.f_B_norm())*h.xi(t,p))

    try: 
        bcs = gr.bcs
    except:
        bcs = [bc_neu, bc_neu]

    d_l, n_l, r_l = find_boundary_condition(bcs, bcleft)
    d_r, n_r, r_r = find_boundary_condition(bcs, bcright)

    if np.all(d_l == d_r) and np.all(n_l == n_r) and np.all(r_l == r_r):
        x = (gr.x - max(gr.x)/2)* h.xi(0,p)/h.xi(t,p)
        xmin = np.min(x)
        xmax = np.max(x)*1.7
    else:
        x = (gr.x)* h.xi(0,p)/h.xi(t,p)
        xmin = 0
        xmax = np.max(x)*1.35
                
        
        
    fig, ax = plt.subplots(2,1, sharex='col')
    # gridspec_kw={'hspace': 0, 'wspace': 0}
    ax[0].plot(x, eden/abs(pot.mat_pars.f_B_norm()) + 1, label="Total")
    ax[0].plot(x, eden_pot/abs(pot.mat_pars.f_B_norm()) + 1, label="Bulk")
    
    ax[0].set_ylabel(r'$e/|f_B|$')
    ax[0].grid()
    ax[0].set_xlim(xmin, xmax)
    ax[0].legend(loc='center right', title=r'Excess over $f_B$')
    ax[0].set_title(r'T={:.2f} mK, p={:.2f} bar: $\sigma/\xi_{{\rm GL}}(T)|f_B(T)| = {:.2f}$'.format(t*h.Tc_mK(p), p, sigma))
    
    norm =  np.sqrt(3)/h.delta_B_norm(t, p) 
    
    for comp in real_comp_list:
        ax[1].plot(x, norm*A[:, comp[0], comp[1]].real, 
                   label=r'${{\rm Re}}(A_{{ {} {} }})$'.format(comp_key[comp[0]], comp_key[comp[1]])
                   )
    for comp in imag_comp_list:
        ax[1].plot(x, norm*A[:, comp[0], comp[1]].imag, 
                   label=r'${{\rm Im}}(A_{{ {} {} }})$'.format(comp_key[comp[0]], comp_key[comp[1]])
                   )
    
    # ax[1].plot(x, norm*A[:, 0, 0].real, label=r'${\rm Re}(A_{xx})$')
    # ax[1].plot(x, norm*A[:, 1, 1].real, label=r'${\rm Re}(A_{yy})$')
    # ax[1].plot(x, norm*A[:, 2, 2].real, label=r'${\rm Re}(A_{zz})$')
    # ax[1].plot(x, norm*A[:, 0, 1].imag, label=r'${\rm Im}(A_{xy})$')
    # ax[1].plot(x, norm*A[:, 1, 0].imag, label=r'${\rm Im}(A_{yx})$')
    
    ax[1].set_xlim(xmin, xmax)
    ax[1].legend(loc='center right')
    ax[1].grid()
    ax[1].set_ylabel(r'$A \sqrt{3}/\Delta_B(T,p)$')
    ax[1].set_xlabel(r'$x/\xi_{\rm GL}(T)$')

    return ax

