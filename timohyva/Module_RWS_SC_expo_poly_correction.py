'''
This script is used for calculating correction of exponent in the RWS SC.

Two temperature scales are possible i.e., PLTS2000 & Greywall1986. There is a 
global symbol table variable Temperature_Scale, which decides the temerature scale.

author: Quang. Zhang (timohyva@github)
'''

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit as cfit

import he3_tools_Vn01 as hn
import Module_SCCO_V03 as SC


#################################################################################
##      >>>>>>>>>>>>>>>>>>   Temperature scale switch  <<<<<<<<<<<<<<<<        ##
#################################################################################

# default value
Temperature_Scale = "PLTS"

def turn_on_Greywall(): 
    '''Function, which could switch Temperature_Scale to "Greywall".

    Using globals() global symbol table dictionary.
    '''
    globals()["Temperature_Scale"] = "Greywall"

def turn_on_PLTS(): 
    '''Function, which could switch Temperature_Scale back to "PLTS" if necessary.

    Using globals() global symbol table dictionary.
    '''
    globals()["Temperature_Scale"] = "PLTS"    

#################################################################################


BO = SC.BETA('betaAndTc')

p_pcp = hn.p_pcp_bar
pressure = np.arange(21.22, 34.2, 0.2) # bar, p_pcp = 21.22
# pressure = np.arange(22., 36., 2) # bar

# weak coupling coefficients of beta_i
cwc = np.array([-1, 2, 2, 2, -2])

nomi = -3.*cwc[0]-cwc[2]+2.*cwc[3]+2*cwc[4]
# print("\n nomi looks like\n", nomi)


##################################################################################
##  >>>>>     model funtions for q_arr fitting with Greywall             <<<<<< ##
##################################################################################

def fit_q_f2_G(x, c0, c1, c2):
    delta_p = x - p_pcp
    return c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f3_G(x, c0, c1, c2, c3):
    delta_p = x - p_pcp
    return c3*(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0
 
def fit_q_f4_G(x, c0, c1, c2, c3, c4):
    delta_p = x - p_pcp
    return c4*(delta_p**4) + c3*(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f5_G(x, c0, c1, c2, c3, c4, c5):
    delta_p = x - p_pcp
    return c5*(delta_p**5) + c4*(delta_p**4) + c3*(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f6_G(x, c0, c1, c2, c3, c4, c5, c6):
    delta_p = x - p_pcp
    return c6*(delta_p**6) + c5*(delta_p**5) + c4*(delta_p**4) + c3*(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0


#############################################################
##        >>>>>>>>>>>>   Greywell 1986  <<<<<<<<<<<<<      ##
#############################################################

q_arr_Greywall = np.array([])

tAB_Greywall = (hn.TAB_poly_Greywall(pressure-hn.p_pcp_bar))/(hn.Tc_poly_Greywall(pressure))

# print("\n tAB_Greywall looks like\n", tAB_Greywall)


for ip in np.arange(0, len(pressure), 1):

    p = pressure[ip]
    
    BO.c1_function(SC.P,SC.c1,p);
    BO.c3_function(SC.P,SC.c3,p);
    BO.c4_function(SC.P,SC.c4,p);
    BO.c5_function(SC.P,SC.c5,p);

    # print("\nnow p is ",p,"\nBO.c1p is ",BO.c1p,"\nBO.c3p is ",BO.c3p,"\nBO.c4p is ",BO.c4p, "\nBO.c5p ",BO.c5p)
    
    dino = 3.*BO.c1p + BO.c3p -2.*BO.c4p -2.*BO.c5p
    # print("\n dino looks like\n", dino)

    tAB_RWS = nomi/dino

    q_arr_Greywall = np.append(q_arr_Greywall, np.log(tAB_RWS/tAB_Greywall[ip]))


# print("\n q_arr with Greywall 1986 scale looks like\n", q_arr_Greywall)

popt2_G, pcov2_G = cfit(fit_q_f2_G, pressure, q_arr_Greywall)
popt3_G, pcov3_G = cfit(fit_q_f3_G, pressure, q_arr_Greywall)
popt4_G, pcov4_G = cfit(fit_q_f4_G, pressure, q_arr_Greywall)
popt5_G, pcov5_G = cfit(fit_q_f5_G, pressure, q_arr_Greywall)
popt6_G, pcov6_G = cfit(fit_q_f6_G, pressure, q_arr_Greywall)

popt_G_list = [popt2_G, popt3_G, popt4_G, popt5_G, popt6_G]

# print("\npopt2_G looks like ",popt2_G,"\npopt3_G looks like ",popt3_G,"\npopt4_G looks like ",popt4_G,"\npopt5_G looks like ",popt5_G,"\npopt6_G looks like ",popt6_G)

# print(" \npopt_G_list looks like ",popt_G_list)

################################################################
################################################################

################################################################
##           >>>>>>>>>>>>   PLTS 2000  <<<<<<<<<<<<<          ##
################################################################


##  >>>>>        model funtions for q_arr fitting with PLTS          <<<<<< ##

def fit_q_f2_PLTS(x, c0, c1, c2): return c2*(x**2) + c1*(x) + c0

def fit_q_f3_PLTS(x, c0, c1, c2, c3): return c3*(x**3) + c2*(x**2) + c1*(x) + c0
 
def fit_q_f4_PLTS(x, c0, c1, c2, c3, c4): return c4*(x**4) + c3*(x**3) + c2*(x**2) + c1*(x) + c0

def fit_q_f5_PLTS(x, c0, c1, c2, c3, c4, c5):
    return c5*(x**5) + c4*(x**4) + c3*(x**3) + c2*(x**2) + c1*(x) + c0

def fit_q_f6_PLTS(x, c0, c1, c2, c3, c4, c5, c6):
    return c6*(x**6) + c5*(x**5) + c4*(x**4) + c3*(x**3) + c2*(x**2) + c1*(x) + c0

################################################################

q_arr_PLTS = np.array([])

tAB_PLTS = (hn.TAB_poly_PLTS(pressure))/(hn.Tc_poly_PLTS(pressure))

# print("\n tAB_PLTS looks like\n", tAB_PLTS)

for ip in np.arange(0, len(pressure), 1):

    p = pressure[ip]
    
    BO.c1_function(SC.P,SC.c1,p);
    BO.c3_function(SC.P,SC.c3,p);
    BO.c4_function(SC.P,SC.c4,p);
    BO.c5_function(SC.P,SC.c5,p);

    # print("\nnow p is ",p,"\nBO.c1p is ",BO.c1p,"\nBO.c3p is ",BO.c3p,"\nBO.c4p is ",BO.c4p, "\nBO.c5p ",BO.c5p)
    
    dino = 3.*BO.c1p + BO.c3p -2.*BO.c4p -2.*BO.c5p
    # print("\n dino looks like\n", dino)

    tAB_RWS = nomi/dino

    q_arr_PLTS = np.append(q_arr_PLTS, np.log(tAB_RWS/tAB_PLTS[ip]))


# print("\n q_arr with PLTS 2000 scale looks like\n", q_arr_PLTS)


popt2_PLTS, pcov2_PLTS = cfit(fit_q_f2_PLTS, pressure, q_arr_PLTS)
popt3_PLTS, pcov3_PLTS = cfit(fit_q_f3_PLTS, pressure, q_arr_PLTS)
popt4_PLTS, pcov4_PLTS = cfit(fit_q_f4_PLTS, pressure, q_arr_PLTS)
popt5_PLTS, pcov5_PLTS = cfit(fit_q_f5_PLTS, pressure, q_arr_PLTS)
popt6_PLTS, pcov6_PLTS = cfit(fit_q_f6_PLTS, pressure, q_arr_PLTS)

popt_PLTS_list = [popt2_PLTS, popt3_PLTS, popt4_PLTS, popt5_PLTS, popt6_PLTS]

# print("\npopt2_PLTS looks like ",popt2_PLTS,"\npopt3 looks like ",popt3_PLTS,"\npopt4 looks like ",popt4_PLTS,"\npopt5 looks like ",popt5_PLTS,"\npopt6 looks like ",popt6_PLTS)

# print(" \npopt_PLTS_list looks like ",popt_PLTS_list)


###################################################################
##   >>>>>>>>>>>>>>>>>>>>     interface      <<<<<<<<<<<<<<<<<<  ##
###################################################################
 
def q(p, n):
    '''Interface for returing the exponent q for fudge coefficient.

    Depending on the value of global variable Temperature_Scale, it return q with 
    PLTS 2000 or Greywall 1986 scales.

    prameter p is pressure in unit of bar;
    in type prameter *n* decides the order of the polynomial; n = 2,3,4,5,6
    '''
    if Temperature_Scale == "PLTS":
        return np.polynomial.polynomial.Polynomial(popt_PLTS_list[n-2])(p) # Polynomial class
    elif Temperature_Scale == "Greywall":
        return np.polynomial.polynomial.Polynomial(popt_G_list[n-2])(p-p_pcp) # Polynomial class
        
        
        
