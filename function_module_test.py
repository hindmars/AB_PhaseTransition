import Module_SCCO_V02 as SC 
import numpy as np

from math import pi
import math

#####################################################################




######################################################################
### build object via SC Object module, and check the intepolatinon ###
######################################################################

BetaObject = SC.BETA('betaAndTc')

######################################################################


######################################################################
#######      Module functions are implemented from here       ########
######################################################################

##

def mEff(p):
   # atomic mass unit, Dalton, kg
   kg = 1;u = 1.66053906660*(10**(-27))*kg 

   # mass of helium3 atom
   m3 = 3.016293*u 

   BetaObject.mStar_function(SC.P,SC.Ms,p); meffective = (BetaObject.ms)*m3
   return meffective
   
# fermi velocity
def vF(p):
   BetaObject.vFermi_function(SC.P,SC.VF,p); vFermi = BetaObject.vf
   return vFermi   

# N(0)
def N0(p):
   J = 1; s = 1;hbar = 1.054571817*(10**(-34))*J*s
   n0 =  ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
#   n0 = J * s * hbar * pi *pi * vF(p) * mEff(p)
   return n0
   

