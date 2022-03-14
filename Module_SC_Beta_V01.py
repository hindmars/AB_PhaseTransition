########################################
# Important notations and descriptions #
########################################


# This script uses the SC BetaObject module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version 0.1 of SC beta parameters code, in which the pico-Joule(pJ 10^-12) unit is used 

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569;

##################################################################

import Module_SCCO_V02 as SC
# import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module
import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math
from matplotlib import ticker, cm


###################################################################
###                  Constants declations                       ###
###################################################################

m = 1;s = 1; J = 1; Kelvin = 1; kg = 1; bar = 1;# Length unit, Time unit, Energy unit, mass unit, pressure unit 
zeta3 = 1.2020569;
# kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1
kb = 1.380649*(10**(-23))*J*(Kelvin**(-1)) # Boltzmann constant J.K^-1
c = 2.99792458*(10**(8))*m*(s**(-1)) # speed of light, m.s^-1

# hbar = 6.582119569*(10**(-16)) #plank constant, eV.s
hbar = 1.054571817*(10**(-34))*J*s # planck constant, J.s
# u = 9.3149410242*(10**(8))*eV*(c**(-2)) # atomic mass unit, Dalton, eV.c^-2
u = 1.66053906660*(10**(-27))*kg # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom


######################################################################
### build object via SC Object module, and check the intepolatinon ###
######################################################################

BetaObject = SC.BETA('betaAndTc')

######################################################################

def calculate_SC_beta_gaps_etc(p,T):

   BetaObject.c1_function(SC.P,SC.c1,p);c1p = BetaObject.c1p
   BetaObject.c2_function(SC.P,SC.c2,p);c2p = BetaObject.c2p
   BetaObject.c3_function(SC.P,SC.c3,p);c3p = BetaObject.c3p
   BetaObject.c4_function(SC.P,SC.c4,p);c4p = BetaObject.c4p
   BetaObject.c5_function(SC.P,SC.c5,p);c5p = BetaObject.c5p
   BetaObject.tc_function(SC.P,SC.Tc,p);Tcp = (BetaObject.tcp)*(10**(-3)) # turn mK to Kelvin 10^(-3)
   BetaObject.mStar_function(SC.P,SC.Ms,p); mEffective = (BetaObject.ms)*m3
   BetaObject.vFermi_function(SC.P,SC.VF,p); vFermi = BetaObject.vf
   BetaObject.xi0_function(SC.P,SC.XI0,p); xi0p = (BetaObject.xi0)*(10**(-9)) # turn nm to m 10^(-9)

   N0 = ((mEffective**(2))*vFermi)/((2*pi*pi)*(hbar**(3))) # energy density of Fermi surface

   t = T/Tcp
    

   print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,' ,xi0p is ', xi0p, '\n\n')

   print('\npressure is, ',p,' effective mass is, ', mEffective, ' Fermi velocity is,', vFermi, ' N(0) is ',N0,'\n\n')
    
   
   print('temperatureis:, ',t)

   xiGLWeakCoupling = math.sqrt((7.0*zeta3)/20.0)*xi0p # weak coupling GL coherent length in PRB 101 024517

   if t >= 1:

      print(" bro, we just got temperature at Tc, you get nothing. **__**, see all stuffs as np.nan ")

      alpha = np.nan
      beta1 = np.nan
      beta2 = np.nan
      beta3 = np.nan
      beta4 = np.nan
      beta5 = np.nan
      betaA = np.nan
      betaB = np.nan
      DeltaA = np.nan
      DeltaB = np.nan
      phi0 = np.nan
      xitp =np.nan
      
                           
   else:

      alpha = (1/3)*N0*(t - 1)
       
      beta_wc1 = -((7*N0*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp));
      beta_wc2 = -2*beta_wc1
      beta_wc3 = -2*beta_wc1
      beta_wc4 = -2*beta_wc1
      beta_wc5 = 2*beta_wc1

      beta_sc1 = c1p*abs(beta_wc1)
      beta_sc2 = c2p*abs(beta_wc1)
      beta_sc3 = c3p*abs(beta_wc1)
      beta_sc4 = c4p*abs(beta_wc1)
      beta_sc5 = c5p*abs(beta_wc1)

      beta1 = beta_wc1 + t*beta_sc1
      beta2 = beta_wc2 + t*beta_sc2
      beta3 = beta_wc3 + t*beta_sc3
      beta4 = beta_wc4 + t*beta_sc4
      beta5 = beta_wc5 + t*beta_sc5 

      xitp = xiGLWeakCoupling/math.sqrt(1-t) # temperature dependent coherent length  
       
      print(' Temperature dependent GL xi is ', xitp)
            
#########################################################
########### betaA, betaB, DeltaA, DeltaB, phi0 ##########
#########################################################

      betaA = beta2 + beta4 + beta5 # betaA

      betaB = beta1 + beta2 + (1/3)*(beta3 + beta4 + beta5) # betaB

      DeltaA = math.sqrt((-alpha)/(2*betaA)) 

      DeltaB = math.sqrt((-alpha)/(2*betaB))

      phi0 = math.sqrt(DeltaA*DeltaA - (math.sqrt(2/3))*DeltaA*DeltaB + DeltaB*DeltaB)
        
   return alpha, beta1, beta2, beta3, beta4, beta5, betaA, betaB, DeltaA, DeltaB, phi0, Tcp, xitp, xi0p, t
           
          






