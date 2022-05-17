############################################################
####        Important notations and descriptions       #####
############################################################


# This script is used for calculating and plotting the ratio f(\phi_{-})/f_{\phi_{+}} of free energy
# deduced from original GL free energy by setting \phi \times D_{\alpha i} = A_{B} - A_{A}. Here A_{A} and A_{B}
# are particular order parameters of A and B phases. The explicit expression see the Mathematica note book.
# The resulted contour plot and density plot show \bar{\lambda} > 0.93 in whole p-T region with H = 0.

# This script uses SC_Beta module to caltulate the M^{2}, \delta and \lambda parameter in the p-T plane,
# then plugging them into the expressions of f(\phi_{-}) and f(\phi_{+}) to get the ratio f(\phi_{-})/f_{\phi_{+}}

# A new version of this code with Lotynk't experimentla parameters will be developed soon based on this code

#log time: 10th. January. 2022


###########################################################################################################
##############         Significant Developments & Processes of WP2 project         ########################
###########################################################################################################

# Two new modules are created based on the primary SC_CorrectionObject module.
# One is SC_Beta module and other one is plot_TAB_line module.

# The basic idea is using A. J. Leggtt's simply considerations when f_A - f_B is small and
# R_c / \xi, R_c/t are both very tiny, and close to T_AB
# For details, checking his Journal of Low Temperature vol.87 571 1992 and Yip & Leggtt review
# in Helium Three 1990.

# the critical radius are plotted as countours, and compared with the experiment data from Lotynk's experiment.
# the result is amazing, and suggests R_c ~ 400 \xi

# Mark suggested another evaluation about R_c in 19th Nov 2021, see the email in details
# the basic idea is R_c ~ ((f_A \xi_GL)/ \Deltaf_AB)/\xi_GL = f_A/ \Deltaf_AB
# temperature region from 0.0 mk to 2.4 mK

# the pressure region is 21 - 34 bar.

# This script uses the SC_Beta module and plot_TAB module, both of them use SC_CorrectionObject module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version.2 of ebergy difference code, in which the pico-Joule(pJ 10^-12) unit is used 

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569;

####################################################################
####                Modules & Constants                        #####                   
####################################################################
##
#'''
#
import Module_SCC_bB_V00 as SCCB # SC_betaB module, dimensionless parameters & gaps
import numpy as np

from math import pi
import math
import csv

#####################################################################
####                  Constants declations                       ####
#####################################################################
##
#'''
#  Those constants are crucial when user want to call physical qulities,
#  such as density of state N(0), symmetry breaking temperature Tc or
#  temperature-dependent GL coherent length.

#'''
#####################################################################

m = 1;s = 1; J = 1; Kelvin = 1; kg =1; bar = 1;# Length unit, Time unit, Energy unit, mass unit, pressure unit 

kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1
c = 2.99792458*(10**(8)) # speed of light, m.s^-1

hbar = 1.054571817*(10**(-34)) # planck constant, J.s

u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom

zeta3 = 1.2020569;

######################################################################
#######      Module functions are implemented from here       ########
######################################################################

##
#'''
#   dimensionless M2, \delta and \lambda coefficients of fphi
#   also \bar{\lambda}
#'''

# dimensionless M2, with unit N(0)
def M2_bar(p, T):
   return ((SCCB.alpha_bar(p, T)
           *((SCCB.alpha_bar(p, T)*(14.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 7.*SCCB.beta3_bar(p, T) + 3.*SCCB.beta4_bar(p, T) + SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T))  
             +8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))
           /(-2.*SCCB.DeltaA2_bar(p, T) + 2.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p,T))*math.sqrt(SCCB.DeltaB2_bar(p, T)) - 2.*SCCB.DeltaB2_bar(p, T)))

# dimensionless \delta, with unit N0 * (Kb * Tc)^(-1)          
def delta3_bar(p, T):          
   return ((3.*math.sqrt(2.)*SCCB.alpha_bar(p, T)
           *((SCCB.alpha_bar(p, T)*(8.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 5.*SCCB.beta3_bar(p, T) + 7.*SCCB.beta4_bar(p, T) + 5.*SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T))  
             +8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))
           /(2.*SCCB.DeltaA2_bar(p, T) - 2.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T)) + 2.*SCCB.DeltaB2_bar(p, T))**1.5)

# domensionless \lambda, with unit N0 * (Kb * Tc)^(-2)
def lambda4_bar(p, T):
   return ((4.*SCCB.alpha_bar(p, T)
           *((SCCB.alpha_bar(p, T)*(5.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 4.*SCCB.beta3_bar(p, T) + 9.*SCCB.beta4_bar(p, T) + 7.*SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T))  
             +8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))
           /(-2.*SCCB.DeltaA2_bar(p, T) + 2.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T)) - 2*SCCB.DeltaB2_bar(p, T))**2) 

# \bar{\lambda}
#'''
# \bar{\lambda} = (9/2) * \lambda * (M2/\delta^{2})
#'''

def lambda_bar(p, T):
   nominator1 = ((SCCB.alpha_bar(p, T)*(14.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 7.*SCCB.beta3_bar(p, T) + 3.*SCCB.beta4_bar(p, T) + SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T)) 
                  +(8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))

   nominator2 = ((SCCB.alpha_bar(p, T)*(5.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 4.*SCCB.beta3_bar(p, T) + 9.*SCCB.beta4_bar(p, T) + 7.*SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T))  
                  +(8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))

   dinominator = ((SCCB.alpha_bar(p, T)*(8.*SCCB.beta1_bar(p, T) + 14.*SCCB.beta2_bar(p, T) + 5.*SCCB.beta3_bar(p, T) + 7.*SCCB.beta4_bar(p, T) + 5.*SCCB.beta5_bar(p, T)))/(3.*SCCB.betaA_bar(p, T)*SCCB.betaB_bar(p, T)) 
                  +(8.*math.sqrt(2./3.)*math.sqrt(SCCB.DeltaA2_bar(p, T))*math.sqrt(SCCB.DeltaB2_bar(p, T))))**2

   return ((nominator1 * nominator2)/dinominator)
   

#######################################################################
####            Defining f(\phi) and f(\phi_{n/p})                 ####
#######################################################################

##
#'''
#  dimensionless f(\phi), \phi should be offered as dimensionless quality
#'''
def fphi_bar(p, T, phi_bar):
   # first term, in unit of N0 * (Kb*Tc)^(2), \phi = \bar{\phi} * (Kb * Tc) 
   f2 = (1./2.) * M2_bar(p, T) * (phi_bar**2)

   # second term, in unit of N0 * (Kb*Tc)^(-1) * (Kb*Tc)^(3)
   f3 = (1./3.) * delta3_bar(p, T) * (phi_bar**3)

   # third term, in unit of N0 * (Kb*Tc)^(-2) * (Kb*Tc)^(4)
   f4 = (1./4.) * lambda4_bar(p, T) * (phi_bar**4)

   return (f2 - f3 + f4)


##
#'''
#  dimensionless f(\phi_n/p), in unit of N0 * (Kb * Tc)^(2)
#'''
def fphi_ebar(p, T):

    f1 = (4./3.) * lambda_bar(p, T)
    f2 = (1 - (8./9.) * lambda_bar(p, T))**(3./2.)
    f3 = (8./27.) * ((lambda_bar(p, T))**2)

    preCoef = ((delta3_bar(p, T))**4)/(24.*(lambda4_bar(p, T))**3)

    fphi_n = preCoef * (-1. + f1 + f2 - f3)
    fphi_p = preCoef * (-1. + f1 - f2 - f3)

    return (fphi_n, fphi_p)


#########################################################################
####                  extrame points \phi_p \phi_n                  #####
#########################################################################

##
#'''
#  dimensionless \phi_n, \phi_p, in unit of Kb*Tc
#'''

def phi_ebar(p, T):
   preCoef = delta3_bar(p, T)/(2.*lambda4_bar(p, T))
    
   phi_n = preCoef * (1. - math.sqrt(1. - (8./9.) * lambda_bar(p, T)))
   phi_p = preCoef * (1. + math.sqrt(1. - (8./9.) * lambda_bar(p, T)))

   return (phi_n, phi_p)






