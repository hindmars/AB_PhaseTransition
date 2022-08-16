#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:20:44 2021

@author: hindmars
"""

"""
Simple routines for finding the critical bubble in a quartic potential. 
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
import matplotlib.pyplot as plt
from scipy.optimize import newton_krylov
import scipy.optimize.nonlin

def get_psi_m(lam_bar):
    return (3 - np.sqrt(9-8*lam_bar))/(2*lam_bar)

def get_psi_b(lam_bar):
    return (3 + np.sqrt(9-8*lam_bar))/(2*lam_bar)


class quartic_paramset:
    """
    Class for the parameters of a quartic potential, including temperature dependence.
    Uses: Enqvist et al Phys Rev D 45, 3415 (1992). Potential is
    $$
    V(\phi,T) = \frac{1}{2}D(T^2 - T_0^2)\phi^2 - \frac{1}{3}ET\phi^3 + \frac{1}{4} \lambda .
    $$
    Note difference in notation for constants $D$, $E$. 

    Method quartic_potential_init returns a tuple suitable for use in creating a quartic_potential object.
    
    A quartic_paramset object can also be directly passed to quartic_potential, along with a temperature.
    """
    def __init__(self, T0, D, E, lam):
        self.T0 = T0
        self.D = D
        self.E = E
        self.lam = lam
#        self.mu0 = mu0 - for future implementation, constant cubic term
        
        self.disc = E**2/(D*lam)
        self.Tc = T0/( 1 - (2/9)*self.disc)

    def m2(self, T):
        """
        Quadratic coefficient
        """
        return  self.D*(T**2 - self.T0**2)
    
    def mu(self, T):
        """
        Cubic coefficient
        """
        return self.E*T
    
    def lam_bar(self, T):
        """
        Useful parameter combination: always between 0 (no barrier) 
        and 1 (critical temperature)
        """
        return (9/2) * self.lam * self.m2(T)/self.mu(T)**2    
        
    def phi_b(self, T):
        """
        Order parameter in broken phase. Returns nan if only one minimum.
        """
        if self.lam_bar(T) <= 9/8:
            val = 0.5*(self.mu(T)/self.lam) * (1 + np.sqrt(1 - (8/9)*self.lam_bar(T)))
        else:
            val = np.nan
            
        return val
            
    def phi_m(self, T):
        """
        Order parameter at maximum of potential. . Returns nan if only one minimum.
        """
        if self.lam_bar(T) <= 9/8:
            val = 0.5*(self.mu(T)/self.lam) * (1 - np.sqrt(1 - (8/9)*self.lam_bar(T)))
        else:
            val = np.nan
            
        return val

    def quartic_potential_init(self, T):
        """
        Returns parameter tuple for use by quartic_potential class.
        """
        return self.phi_m(T), self.lam, self.phi_b(T)


class quartic_potential:
    """ 
    Quartic potential class. Creates a quartic potential object with a tuple of   
    three necessary parameters (see v function) 
    m2 - quadratic parameter (mass sqaured)
    mu - cubic parameter
    lam - quartic paramater
    where
        V(phi) = V0 + 0.5 * m2 * phi**2 - (1/3)*mu*phi**3 + 0.25*lam*phi**4
    (minumum value of potential is always zero, always an extremum at phi = 0)
    Stores as class vaiables
        phi_b - broken phase field value, default 1
        phi_m - unstable maximum phi value

    These values can be supplied by 
    - quartic_paramset object, in which case temperature T must be supplied.
    - lambda_bar (float) only (lam = 1, phi_b = 1)
    - phi_m (float) only (lam = 1, phi_b = 1)
    """
    
    def __init__(self, params, T=None, lambda_bar=True):
        if isinstance(params, tuple):
            self.m2 = params[0]
            self.mu = params[1]
            self.lam = params[2]
            self.lam_bar = (9/2) * self.m2/self.mu**2
            self.phi_m = (self.mu/(2*self.lam))*(1 - np.sqrt(1 - (8/9)*self.lam_bar))
            self.phi_b = (self.mu/(2*self.lam))*(1 + np.sqrt(1 - (8/9)*self.lam_bar))
        elif isinstance(params, quartic_paramset):
            self.m2 = params.m2(T)
            self.mu = params.mu(T)
            self.lam_bar = (9/2) * self.m2/self.mu**2
            self.phi_m = params.phi_m(T)
            self.lam = params.lam
            self.phi_b = params.phi_b(T)
        elif isinstance(params, float):
            if lambda_bar:
                self.phi_m = get_psi_m(params)/get_psi_b(params)
                self.lam = 1.0
                self.phi_b = 1.0
            else:
                self.phi_m = params
                self.lam = 1.0
                self.phi_b = 1.0
                                
            self.m2 = self.lam * self.phi_m * self.phi_b
            self.mu = self.lam * (self.phi_b + self.phi_m) 
            self.lam_bar = (9/2) * self.m2/self.mu**2

        else:
            raise ValueError("quartic_potential.__init__(): input type not recognised")

        self.T = T
        self.lam_phi_b_4 = self.lam * self.phi_b**4
        self.x_m = self.phi_m/self.phi_b
        self.v0 = (1/12) * self.lam * self.phi_b**3 * (self.phi_b - 2*self.phi_m)


    def v(self, phi):
        # x = phi/self.phi_b
        # return 0.25*self.lam_phi_b_4*(x**2*(x-1)**2  + (1/3)*(1 - 2*self.x_m)*(1 - 3*x**2 + 2*x**3))
        return self.v0 + 0.5 * self.m2 * phi**2 - (1/3)*self.mu*phi**3 + 0.25*self.lam*phi**4
    
    def v_prime(self, phi):
        # x = phi/self.phi_b
        # return self.lam_phi_b_4*x*(x - self.x_m)*(x-1)
        return self.m2 * phi - self.mu*phi**2 + self.lam*phi**3
    
    def v_prime_prime(self,phi):
        # x = phi/self.phi_b
        # return self.lam_phi_b_4*((x - self.x_m)*(x-1) + x*(x-1) + x*(x - self.x_m))
        return self.m2 - 2*self.mu*phi + 3*self.lam*phi**2
    
    def plot(self):
        phi = np.linspace(0,1.2*self.phi_b,100)
        fig, ax = plt.subplots()
        ax.plot(phi/self.phi_b, self.v(phi))
        ax.set_xlabel(r'$\phi/\phi_{\rm b}$')
        ax.set_ylabel(r'$V(\phi)$')
        ax.grid()
        return ax
        

    
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
        self.surface_energy = np.sqrt(pot.lam/2) * pot.phi_b**3/6
        self.energy_diff = pot.v(0)
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
                
        else:
            self.r_bub = np.inf
            self.e_bub = np.inf


def field_eqn(phi, pot, gr):
    """
    Scalar field equation in 1 radial dimension, with trivial gradient term.
    """
    
    phi_plus = np.roll(phi,-1)
    phi_minus = np.roll(phi,+1)

    # Neumann boundary conditions - derivative vanishes at r=0, max(r)    
    phi_plus[-1] = phi[-2]
    phi_minus[0] = phi[1]
    
    second = (phi_plus + phi_minus - 2*phi)/gr.dx2
    first = (phi_plus  - phi_minus)/(gr.x + gr.x_minus)
    first *= (gr.dim - 1)/gr.dx
    
    return second + first - pot.v_prime(phi)


def initial_condition(gr, pot, bub):
    """
    Initial guess: a kink-antikink pair symmetric around r=0, at +/- r_bub
    """
    m2 = pot.v_prime_prime(pot.phi_b)
    k = np.sqrt(m2/4)

    rb = bub.r_bub
    if rb == np.inf:
        rb = gr.R/2

    phi_init = 0.25*(1 - np.tanh(k*(gr.x - rb)))*(1 + np.tanh(k*(gr.x + rb )))  

    return phi_init


def krylov_bubble(param, T=None, gr_pars=(200,20), dim=3, lambda_bar=True, display=False):
    """
    Apply Krylov solver to find unstable bubble solution. 
    If temperature T not given, assumes param is value of order parameter at maximum, for 
    quartic potential (with lam = 1, phi_b = 1).
    If temperature T is given, param is a quartic_paramset, to be evaluated at temperature T.
    
    Returns: order parameter phi, the potential object, and the grid object.
    """
    
    if T is not None:
        phi_m = param.phi_m(T)
        phi_b = param.phi_b(T)
        pot = quartic_potential(param, T)
    else:
        phi_b = 1.0
        pot = quartic_potential(param, lambda_bar=lambda_bar)
        
    
    # print(type(phi_m), phi_m)
    gr = grid_1d(*gr_pars, dim=dim)
    bub = thin_wall_bubble(pot, dim=dim)
    
    phi_init = initial_condition(gr, pot, bub)
    
    def field_eqn_fix(phi):
        return field_eqn(phi, pot, gr)

    try:
        phi = newton_krylov(field_eqn_fix, phi_init)
    except scipy.optimize.nonlin.NoConvergence as e:
        phi = e.args[0]
        print('No Convergence')

    if display:
        
        plt.figure()
        plt.plot(gr.x, phi_init/phi_b, label='Thin wall') 
        plt.plot(gr.x, phi/phi_b, label='Final' )
        plt.xlabel(r'$r$')
        plt.ylabel(r'$\phi/\phi_{\rm b}$')
        plt.xlim(0,gr.R)
        plt.grid()
        plt.legend()

    return phi, pot, gr 


def energy(phi, pot, gr):
    phi_plus = np.roll(phi,-1)
    phi_plus[-1] = phi[-2]
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
    
    first = (phi_plus  - phi)/gr.dx
    e_grad = vol_factor * 0.5*np.trapz(first**2, gr.x**(gr.dim) )
    e_pot = vol_factor * np.trapz(pot.v(phi)-pot.v(0), gr.x**(gr.dim))
    
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


