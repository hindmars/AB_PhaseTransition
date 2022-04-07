
######################################################
'''
this script is used for calculating correction of exponent in the RWS SC

author: Quang (timohyva@github)
'''

import numpy as np
import math
import matplotlib as plt

import he3_tools_Vn01 as h
import Module_SCCO_V02 as SC

BO = SC.BETA('betaAndTc')

pressure = np.arange(21.22, 34.2, 0.2) # bar

lambda_arr = np.array([])

tAB = (h.TAB_poly_Greywall(pressure-h.p_pcp_bar))/(h.Tc_poly_Greywall(pressure))

print("\n tAB looks like\n", tAB)

# weak coupling coefficients
cwc = np.array([-1, 2, 2, 2, -2])

nomi = -3.*cwc[0]-cwc[2]+2.*cwc[3]+2*cwc[4]

print("\n nomi looks like\n", nomi)

for ip in np.arange(0, len(pressure), 1):

    p = pressure[ip]
    
    BO.c1_function(SC.P,SC.c1,p);
    BO.c3_function(SC.P,SC.c2,p);
    BO.c4_function(SC.P,SC.c3,p);
    BO.c5_function(SC.P,SC.c3,p);
    
    dino = 3.*BO.c1p + BO.c3p -2.*BO.c4p -2.*BO.c5p
    print("\n dino looks like\n", dino)

    lambda_arr = np.append(lambda_arr, (math.log(nomi/dino))/(math.log(tAB[ip])))

print("\n lambda_arr looks like\n", lambda_arr)
