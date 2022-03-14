####################################################################
################ Important notations and descriptions ##############
####################################################################
#
#                   >>>>  READ ME ! Please !  <<<<
#

# this is version 0.0 of module of Strong Coupling (SC) correction. which provides easy-using dimensionless material parameters, gaps.

# It provides SC corrected Ginzburg-Landau coefficients, gap energies and other physical qulities
# with pattern 
     
#    " dimentionless expressions * dimensional expressions ".

# In this way, user can get dimensionless form of physical qualities by using object attribute pattern.
#
# It is based on another module called Module_SC_Beta_V0*.py, which can return dimensionless parameters with very detailed way.
# When the dimenssional expressions of these qualities are intrested, user can import that module as well.

# This module will be suitable for most of regular calculations, in which user won't involve huge exponent from N(0).
# And all paramters and qualities can be called as object attributes after importing.
   


# For example, version 0.1 gives out \bar{\alpha} and \bar{\beta_{1}} as 
#'''
#  alpha_bar (in unit of  N(0)],

#  beta_bar  (in unit of N(0) (kb * Tc)^(-2))

#  Delta2_bar (\Delta^{2} in unit of (kb*Tc)^(2)) 

#'''
# pnoa(b), N(0), Tc and betatd will be offered separately when they are called with " module.name " form i.e., SSC.beta_td etc.

####################################################################

# This script uses the SC_Beta_V0*.py module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569

#####################################################################

import Module_SC_Beta_V03 as SCtd 
import numpy as np

from math import pi
import math

#####################################################################
####                  Constants declations                       ####
#####################################################################
####
#'''
#  Those constants are crucial when user want to call physical qulities,
#  such as density of state N(0), symmetry breaking temperature Tc or
#  temperature-dependent GL coherent length.

#'''
#####################################################################

# Length unit, Time unit, Energy unit, mass unit, pressure unit 
m = 1;s = 1; J = 1; Kelvin = 1; kg = 1; bar = 1;

# \zeta(3)
zeta3 = 1.2020569;

# Boltzmann constant J.K^-1
kb = 1.380649*(10**(-23))*J*(Kelvin**(-1))

# speed of light, m.s^-1
c = 2.99792458*(10**(8))*m*(s**(-1)) 

# planck constant, J.s
hbar = 1.054571817*(10**(-34))*J*s 

# atomic mass unit, Dalton, kg
u = 1.66053906660*(10**(-27))*kg 

# mass of helium3 atom
m3 = 3.016293*u 

######################################################################


######################################################################
#######      Module functions are implemented from here       ########
######################################################################

##
#'''
# Tcp, pnoa, pnob, N(0) 
#'''
# Tc in mK
def Tcp(p): return SCtd.Tcp(p)
   
# effective mass in kg
def mEff(p): return SCtd.mEff(p)
   
# fermi velocity
def vF(p): return SCtd.vF(p)
   
# xi0p in nanometer
def xi0p(p): return SCtd.xi0p(p)
   
# N(0)
def N0(p): return SCtd.N0(p)
                                                          

######################################################################
##
# \bar{\alpha}, \bar{\beta}_{i} in unit of N0, [N(0) (kb * Tc)^(-2)] respctive
#                                 

# dimensionless alpha i,e., \bar{\alpha} or alpha_bar
def alpha_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_alpha*SCtd.alpha_td(p, T)
                                                                                      
# dimensionless beta_1 i.e., \bar{{\beta}}_1 or beta1_bar
def beta1_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.beta1_td(p, T)

# dimensionless beta_2 i.e., \bar{\beta}_2  or beta2_bar
def beta2_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.beta2_td(p, T)

# dimensionless beta_3 i.e., \bar{\beta}_3 or beta3_bar
def beta3_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.beta3_td(p, T)

# dimensionless beta_4 i.e., \bar{\beta}_4 or beta4_bar
def beta4_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.beta4_td(p, T)

# dimensionless beta_5 i.e., \bar{\beta}_5 or beta5_bar
def beta5_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.beta5_td(p, T)


########################################################################   
##
#  \bar{\beta}_{A}, \bar{\beta}_{B} in unit of [N(0) (kb * Tc)^(-2)]
#                                

# \beta_A  i.e., \bar{\beta}_A or betaA_bar                                
def betaA_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.betaA_td(p, T)

# \beta_B  i.e., \bar{\beta}_B or betaB_bar                                                               
def betaB_bar(p, T):
   (c_alpha, c_betai) = SCtd.pno(); return c_betai*SCtd.betaB_td(p, T)


########################################################################
##
#  \Delta_{A}^{2}, \Delta_{B}^{2} in unit of (kb * Tc)^2
#

# \Delta_{A}^{2} i.e., \bar{\Delta_{A}^{2}} or DeltaA2_bar
def DeltaA2_bar(p, T): return -(alpha_bar(p, T)/(2.0*betaA_bar(p, T)))

# \Delta_{B}^{2} i.e., \bar{\Delta_{B}^{2}} or DeltaB2_bar
def DeltaB2_bar(p, T): return -(alpha_bar(p, T)/(2.0*betaB_bar(p, T)))
