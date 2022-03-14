####################################################################
################ Important notations and descriptions ##############
####################################################################
#
#                   >>>>  READ ME ! Please !  <<<<
#

# this is version 0.2 of module " Strong Coupling (SC) correction ".

# It provides SC corrected Ginzburg-Landau coefficients, gap energies and other physical qulities
# with tamplate 
     
#    " pno * dimentionless expressions * dimensional expressions ".

# In this way, user can get dimensionless form of physical qualities by using object attribute pattern.
# When the dimenssional expressions of these qualities are intrested, user can easily combine the pre-index (pno)
# and the dimensional expressions to the dimensionless expression to get the desired qualities.

# pno and dimensional expressions are also be provided by this module, and can be called as object attributes.
   


# For example, version 0.2 gives out \alpha and \beta_{1} as 
#'''
#   pnoa * alpha_td * N(0),

#   pnob * beta_td * N(0) (kb * Tc)^(-2).

#'''
# pnoa(b), N(0), Tc and betatd will be offered separately when they are called with " module.name " form i.e., SSC.beta_td etc.

####################################################################

# This script uses the SC BetaObject module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569

#####################################################################

import Module_SCCO_V02 as SC 
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
### build object via SC Object module, and check the intepolatinon ###
######################################################################

BetaObject = SC.BETA('betaAndTc')

######################################################################


######################################################################
#######      Module functions are implemented from here       ########
######################################################################

##
#'''
# Tcp, pnoa, pnob, N(0) 
#'''
# Tc in mK
def Tcp(p):
   Kelvin = 1.; mK = 10**(-3)
   BetaObject.tc_function(SC.P,SC.Tc,p);

   return (BetaObject.tcp)*mK

# effective mass
def mEff(p):
   # atomic mass unit, Dalton, kg
   kg = 1;u = 1.66053906660*(10**(-27))*kg 

   # mass of helium3 atom
   m3 = 3.016293*u 

   BetaObject.mStar_function(SC.P,SC.Ms,p); 

   return (BetaObject.ms)*m3

# fermi velocity
def vF(p):
   BetaObject.vFermi_function(SC.P,SC.VF,p); return BetaObject.vf

# xi0p
def xi0p(p):
   # nano meter
   m = 1; nm = (10**(-9))*m 
   
   BetaObject.xi0_function(SC.P,SC.XI0,p); return (BetaObject.xi0)*nm

# N(0)
def N0(p):
   J = 1; s = 1;hbar = 1.054571817*(10**(-34))*J*s

   return ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
   
# pnoa and pnob tuple
def pno(): zeta3 = 1.2020569; return (1./3., (7.0*zeta3)/(240.0*pi*pi))

   
                                                         

######################################################################
##
# \tilde{\alpha}, \tilde{\beta}_{i}
#                                 

# dimensionless alpha i,e., \tilde{\alpha}
def alpha_td(p, T): return (T/Tcp(p) - 1.) 

                                                                                      
# dimensionless beta_1 i.e., \tilde{\beta}_1
def beta1_td(p, T):
   BetaObject.c1_function(SC.P,SC.c1,p); return -1. + (T/Tcp(p))*BetaObject.c1p

# dimensionless beta_2 i.e., \tilde{\beta}_2 
def beta2_td(p, T):
    BetaObject.c2_function(SC.P,SC.c2,p); return 2. + (T/Tcp(p))*BetaObject.c2p

# dimensionless beta_3 i.e., \tilde{\beta}_3
def beta3_td(p, T):
   BetaObject.c3_function(SC.P,SC.c3,p); return  2. + (T/Tcp(p))*BetaObject.c3p

# dimensionless beta_4 i.e., \tilde{\beta}_4
def beta4_td(p, T):
   BetaObject.c4_function(SC.P,SC.c4,p); return 2. + (T/Tcp(p))*BetaObject.c4p

# dimensionless beta_5 i.e., \tilde{\beta}_5
def beta5_td(p, T):
   BetaObject.c5_function(SC.P,SC.c5,p); return -2. + (T/Tcp(p))*BetaObject.c5p


########################################################################   
##
#  \tilde{\beta}_{A}, \tilde{\beta}_{B}
#                                

# \beta_A                                 
def betaA_td(p, T): return beta2_td(p, T) + beta4_td(p, T) + beta5_td(p, T)

# \beta_B                                 
def betaB_td(p, T):
    return (beta1_td(p, T) + beta2_td(p, T)
            + (1./3.)*(beta3_td(p, T) + beta4_td(p, T) + beta5_td(p, T)))
