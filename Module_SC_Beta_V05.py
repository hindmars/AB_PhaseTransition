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

import Module_SCCO_V03 as SC 
import numpy as np

from math import pi
import math

import Module_RWS_SC_expo_poly_correction_Greywall as SC_CoEpoly_G

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
zeta3 = 1.2020569031595942854

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
###   string for switch between expo-correccted SC & original SC   ###
######################################################################

SC_Correction_Switch = "ON"

# SC_Correction_Switch = "OFF" 

######################################################################

######################################################################
##'''
#

# WimanSauls15_SC = "YES"
WimanSauls15_SC = "NO"

######################################################################

######################################################################
##'''        fudge coefficients array
#

fc_arr = np.array([1., 1., 1., 1., 1.])

#####################################################################


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

# xiGL(p,T) from JWS PRB, there is wired " 20 " 
def xiGL_JWS(p, T):
   # p,T dependent GL coherent length, nanometer
    zeta3 = 1.2020569031595942854; xiGL = xi0p(p)*(math.sqrt((7.*zeta3)/20.));
    return xiGL/(math.sqrt(1.-T/Tcp(p)))
 
# xiGL(p, T) from Sauls Mithushima paper and Osheoff Cross 1977 PRL about mesument of surface tension
def xiGL_OC(p, T):
   # p,T dependent GL coherent length, nanometer
    zeta3 = 1.2020569031595942854; xiGL = xi0p(p)*(math.sqrt((7.*zeta3)/12.));
    return xiGL/(math.sqrt(1.-T/Tcp(p)))
 
# N(0)
def N0(p):
   J = 1; s = 1;hbar = 1.054571817*(10**(-34))*J*s

   return ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))

# K1, K2, K3, and K for gradient terms
def K(p):
   zeta3 = 1.2020569031595942854;
   return ((7.*zeta3)/60.)*N0(p)*xi0p(p)*xi0p(p)   
   
# pnoa and pnob tuple
def pno(): zeta3 = 1.2020569; return (1./3., (7.0*zeta3)/(240.0*pi*pi))

   
                                                         

######################################################################
##
# \tilde{\alpha}, \tilde{\beta}_{i}
#                                 

# dimensionless alpha i,e., \tilde{\alpha}
def alpha_td(p, T): return (T/Tcp(p) - 1.) 

                                                                                      
# dimensionless beta_1 i.e., \tilde{\beta}_1
def beta1_td(p, T, fudge_c=fc_arr, key=SC_Correction_Switch, scOption_WS=WimanSauls15_SC):
   BetaObject.c1_function(SC.P,SC.c1,p);

   # print("\n key looks like ", key)
   
   if (key == "ON") and (scOption_WS == "NO"):
     # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
     # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
     # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
     # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
     q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4)  

     # print("\n q looks like ", q)

     if q < 0:
       return -1. + (T/Tcp(p))*np.exp(q)*BetaObject.c1p
     elif q >= 0:
       # return -1. + (T/Tcp(p))*BetaObject.c1p
        return -1. + (T/Tcp(p))*np.exp(q)*BetaObject.c1p
    
   elif (key  == "OFF")  and (scOption_WS == "NO"):
      # return -1. + (T/Tcp(p))*BetaObject.c1p
      # print(" fudge_c looks like\n ", fudge_c,"\n")
      return -1. + (T/Tcp(p))*(fudge_c[0])*BetaObject.c1p

   if scOption_WS == "YES":
      return -1. + (T/Tcp(p))*SC.b1sc_poly(p)

# dimensionless beta_2 i.e., \tilde{\beta}_2 
def beta2_td(p, T, fudge_c=fc_arr, key=SC_Correction_Switch, scOption_WS=WimanSauls15_SC):
    BetaObject.c2_function(SC.P,SC.c2,p);

    if (key == "ON") and (scOption_WS == "NO"):
      # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
      # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
      # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
      # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
      q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4)

      if q < 0:
        return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c2p
      elif q >= 0:
         # return 2. + (T/Tcp(p))*BetaObject.c2p
         return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c2p
    
    elif (key  == "OFF") and (scOption_WS == "NO"):
       # return 2. + (T/Tcp(p))*BetaObject.c2p
       # print(" fudge_c looks like\n ", fudge_c,"\n")
       return 2. + (T/Tcp(p))*(fudge_c[1])*BetaObject.c2p

    if scOption_WS == "YES":
      return 2. + (T/Tcp(p))*SC.b2sc_poly(p)


    # return 2. + (T/Tcp(p))*BetaObject.c2p

# dimensionless beta_3 i.e., \tilde{\beta}_3
def beta3_td(p, T, fudge_c=fc_arr, key=SC_Correction_Switch, scOption_WS=WimanSauls15_SC):
   BetaObject.c3_function(SC.P,SC.c3,p);

   if (key == "ON") and (scOption_WS == "NO"):
      # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
      # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
      # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
      # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
      q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4)

      if q < 0:
        return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c3p
      elif q >= 0:
         # return 2. + (T/Tcp(p))*BetaObject.c3p
         return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c3p
    
   elif (key == "OFF") and (scOption_WS == "NO"):
       # return 2. + (T/Tcp(p))*BetaObject.c3p
       # print(" fudge_c looks like\n ", fudge_c,"\n")
       return 2. + (T/Tcp(p))*(fudge_c[2])*BetaObject.c3p
      
   if scOption_WS == "YES":
     return 2. + (T/Tcp(p))*SC.b3sc_poly(p)
   
   # return  2. + (T/Tcp(p))*BetaObject.c3p

# dimensionless beta_4 i.e., \tilde{\beta}_4
def beta4_td(p, T, fudge_c=fc_arr, key=SC_Correction_Switch, scOption_WS=WimanSauls15_SC):
   BetaObject.c4_function(SC.P,SC.c4,p);

   if (key == "ON") and (scOption_WS == "NO"):
      # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
      # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
      # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
      # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
      q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4) 

      if q < 0:
        return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c4p
      elif q >= 0:
         # return 2. + (T/Tcp(p))*BetaObject.c4p
         return 2. + (T/Tcp(p))*np.exp(q)*BetaObject.c4p
    
   elif (key == "OFF") and (scOption_WS == "NO"):
       #return 2. + (T/Tcp(p))*BetaObject.c4p
       #print(" fudge_c looks like\n ", fudge_c,"\n")
       return 2. + (T/Tcp(p))*(fudge_c[3])*BetaObject.c4p

   if scOption_WS == "YES":
     return 2. + (T/Tcp(p))*SC.b4sc_poly(p)

   # return 2. + (T/Tcp(p))*BetaObject.c4p

# dimensionless beta_5 i.e., \tilde{\beta}_5
def beta5_td(p, T, fudge_c=fc_arr, key=SC_Correction_Switch, scOption_WS=WimanSauls15_SC):
   BetaObject.c5_function(SC.P,SC.c5,p);

   if (key == "ON") and (scOption_WS == "NO"):
      # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
      # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
      # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
      # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
      q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4)

      if q < 0:
        return -2. + (T/Tcp(p))*np.exp(q)*BetaObject.c5p
      elif q >= 0:
         # return -2. + (T/Tcp(p))*BetaObject.c5p
         return -2. + (T/Tcp(p))*np.exp(q)*BetaObject.c5p
    
   elif (key == "OFF") and (scOption_WS == "NO"):
       #return -2. + (T/Tcp(p))*BetaObject.c5p
       print(" fudge_c looks like\n ", fudge_c,"\n")
       return -2. + (T/Tcp(p))*(fudge_c[4])*BetaObject.c5p

   if scOption_WS == "YES":
     return -2. + (T/Tcp(p))*SC.b5sc_poly(p)

   # return -2. + (T/Tcp(p))*BetaObject.c5p


########################################################################   
##
#  \tilde{\beta}_{A}, \tilde{\beta}_{B}, fA, fB, \Delta_{A}, \Delta_{B}
#                                

# \beta_A                                 
def betaA_td(p, T, fudge_c=fc_arr): return beta2_td(p, T, fudge_c) + beta4_td(p, T, fudge_c) + beta5_td(p, T, fudge_c)

# \beta_B                                 
def betaB_td(p, T, fudge_c=fc_arr):
    return (beta1_td(p, T, fudge_c) + beta2_td(p, T, fudge_c)
            + (1./3.)*(beta3_td(p, T, fudge_c) + beta4_td(p, T, fudge_c) + beta5_td(p, T, fudge_c)))

# bulk fA, with order parameter (\Delta/sqrt(2))*x_alpha*(x_i + i*y_i)
def fA(p, T):
   '''Bulk GL free energy of uniform A phase.

   order parameter is (\Delta/sqrt(2))*x_alpha*(x_i + i*y_i)'''
   J = 1.; Kelvin = 1.
   kb = 1.380649*(10**(-23))*J*(Kelvin**(-1))

   pno_tuple = pno();pnoa = pno_tuple[0];pnob = pno_tuple[1]

   coef = -(pnoa**2)/(4.*pnob)

   return coef*((kb*Tcp(p))**2)*N0(p)*((alpha_td(p, T)*alpha_td(p, T))/betaA_td(p, T))

# bulk fB, with order parameter (\Delta/sqrt(3))*\delta_{a,i}
def fB(p, T):
   '''Bulk GL free energy of uniform B phase.

   order parameter is (\Delta/sqrt(3))*\delta_{a,i}'''
   J = 1.; Kelvin = 1.
   kb = 1.380649*(10**(-23))*J*(Kelvin**(-1))

   pno_tuple = pno();pnoa = pno_tuple[0];pnob = pno_tuple[1]

   coef = -(pnoa**2)/(4.*pnob)

   return coef*((kb*Tcp(p))**2)*N0(p)*((alpha_td(p, T)*alpha_td(p, T))/betaB_td(p, T))

def GapA(p, T):
   '''\Delta_{A}(p, T) for uniform A phase.'''
   J = 1.; Kelvin = 1.
   kb = 1.380649*(10**(-23))*J*(Kelvin**(-1))

   pno_tuple = pno();pnoa = pno_tuple[0];pnob = pno_tuple[1]

   coef = -(1./2.)*(pnoa/pnob)

   return np.sqrt(coef*((kb*Tcp(p))**2)*(alpha_td(p, T)/betaA_td(p, T)))

def GapB(p, T):
   '''\Delta_{B}(p, T) for uniform B phase.'''
   J = 1.; Kelvin = 1.
   kb = 1.380649*(10**(-23))*J*(Kelvin**(-1))

   pno_tuple = pno();pnoa = pno_tuple[0];pnob = pno_tuple[1]

   coef = -(1./2.)*(pnoa/pnob)

   return np.sqrt(coef*((kb*Tcp(p))**2)*(alpha_td(p, T)/betaB_td(p, T)))

   


########################################################################
##
#  
#


# tAB from RWS19 SC
def tAB_RWS(p):
   BetaObject.c1_function(SC.P,SC.c1,p);
   BetaObject.c3_function(SC.P,SC.c3,p);
   BetaObject.c4_function(SC.P,SC.c4,p);
   BetaObject.c5_function(SC.P,SC.c5,p);

   tab_rws = 1./(3.*BetaObject.c1p + BetaObject.c3p - 2.*BetaObject.c4p - 2.*BetaObject.c5p) 

   if (tab_rws <= 1.) and (tab_rws > 0.):
      return tab_rws
   else:
      return np.nan

   # return tab_rws

def tAB_RWSfarr(p, fudge_c=fc_arr):
   BetaObject.c1_function(SC.P,SC.c1,p);
   BetaObject.c3_function(SC.P,SC.c3,p);
   BetaObject.c4_function(SC.P,SC.c4,p);
   BetaObject.c5_function(SC.P,SC.c5,p);

   tab_rwsfarr = 1./(3.*(fudge_c[0])*BetaObject.c1p + (fudge_c[2])*BetaObject.c3p - 2.*(fudge_c[3])*BetaObject.c4p - 2.*(fudge_c[4])*BetaObject.c5p) 

   if (tab_rwsfarr <= 1.) and (tab_rwsfarr > 0.):
      return tab_rwsfarr
   else:
      return np.nan

def tAB_WS15_interp(p):
   BetaObject.c1_function(SC.P,SC.c1,p);
   BetaObject.c3_function(SC.P,SC.c3,p);
   BetaObject.c4_function(SC.P,SC.c4,p);
   BetaObject.c5_function(SC.P,SC.c5,p);

   tab_rws = 1./(3.*BetaObject.c1p + BetaObject.c3p - 2.*BetaObject.c4p - 2.*BetaObject.c5p) 

   if (tab_rws <= 1.) and (tab_rws > 0.):
      return tab_rws
   else:
      return np.nan   

def tAB_WS15(p):

   # tab_ws = 1./(3.*BetaObject.c1p + BetaObject.c3p - 2.*BetaObject.c4p - 2.*BetaObject.c5p)
   tab_ws = 1./(3.*SC.b1sc_poly(p) + SC.b3sc_poly(p) - 2.*SC.b4sc_poly(p) - 2.*SC.b5sc_poly(p))
   # SC.b5sc_poly(p)

   if (tab_ws <= 1.) and (tab_ws > 0.):
      return tab_ws
   else:
      return np.nan
   
     

# correccted tAB from RWS19 SC, only work when SC_Correction_Switch == "NO"
def tAB_RWSco(p):
   
   if SC_Correction_Switch == "ON":
   #if SC_Correction_Switch == "OFF":
      # q = SC_CoEpoly_G.fit_q_f0x2(p, *SC_CoEpoly_G.popt0x2)
      # q = SC_CoEpoly_G.fit_q_f2(p, *SC_CoEpoly_G.popt2)
      # q = SC_CoEpoly_G.fit_q_f3(p, *SC_CoEpoly_G.popt3)
      # q = SC_CoEpoly_G.fit_q_f6(p, *SC_CoEpoly_G.popt6)
      q = SC_CoEpoly_G.fit_q_f4(p, *SC_CoEpoly_G.popt4)
      # q = 0.9

      tab_rwsco = np.exp(-q)*tAB_RWS(p)
      # tab_rwsco = (q)*tAB_RWS(p)

      if (tab_rwsco <= 1.) and (tab_rwsco > 0.):
         return tab_rwsco
      else:
         return np.nan
      # else:
      #    return np.nan
      
      # if np.isnan(tAB_RWS(p)) == True:
      #    return np.nan  
      # elif: (tab_rwsco <= 1.) and (tab_rwsco > 0.):

      #    tab_rwsco = math.exp(-q)*tAB_RWS(p)
      #    return np.nan

# return TAB_RWS in unit of Kelvin
def TAB_RWS(p):
   # if np.isnan(tAB_RWS(p)) == True :
   #    return np.nan
   # else:
   #    return tAB_RWS(p) * Tcp(p)

   return tAB_RWS(p) * Tcp(p)

# return TAB_RWSco in unit of Kelvin, only work when SC_Correction_Switch == "NO"
def TAB_RWSco(p):
   # if np.isnan(tAB_RWSco(p)) == True:
   #   return np.nan
   # else:
     return tAB_RWSco(p) * Tcp(p)

# return TAB_WS15 in unit of Kelvin
def TAB_WS15(p):
   
   return tAB_WS15(p) * Tcp(p)

# return TAB_WS15 interplation in unit of Kelvin
def TAB_WS15_interp(p):
   
   return tAB_WS15_interp(p) * Tcp(p)

