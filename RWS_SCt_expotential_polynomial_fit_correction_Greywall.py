
######################################################
'''
this script is used for calculating correction of exponent in the RWS SC

author: Quang (timohyva@github)
'''

import numpy as np
import math
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit as cfit

import he3_tools_Vn01 as h
import Module_SCCO_V02 as SC

BO = SC.BETA('betaAndTc')

p_pcp_Greywall = h.p_pcp_bar
pressure = np.arange(21.22, 34.2, 0.2) # bar, p_pcp = 21.22
# pressure = np.arange(22., 36., 2) # bar
p_arr = np.arange(20.0, 34.2, 0.2) # bar, p_pcp = 21.22

q_arr = np.array([])

tAB_Greywall = (h.TAB_poly_Greywall(pressure-h.p_pcp_bar))/(h.Tc_poly_Greywall(pressure))

print("\n tAB looks like\n", tAB_Greywall)

# weak coupling coefficients
cwc = np.array([-1, 2, 2, 2, -2])

nomi = -3.*cwc[0]-cwc[2]+2.*cwc[3]+2*cwc[4]

print("\n nomi looks like\n", nomi)

for ip in np.arange(0, len(pressure), 1):

    p = pressure[ip]
    
    BO.c1_function(SC.P,SC.c1,p);
    BO.c3_function(SC.P,SC.c3,p);
    BO.c4_function(SC.P,SC.c4,p);
    BO.c5_function(SC.P,SC.c5,p);

    print("\nnow p is ",p,"\nBO.c1p is ",BO.c1p,"\nBO.c3p is ",BO.c3p,"\nBO.c4p is ",BO.c4p, "\nBO.c5p ",BO.c5p)
    
    dino = 3.*BO.c1p + BO.c3p -2.*BO.c4p -2.*BO.c5p
    print("\n dino looks like\n", dino)

    tAB_RWS = nomi/dino

    q_arr = np.append(q_arr, np.log(tAB_RWS/tAB_Greywall[ip]))

    # lambda_arr = np.append(lambda_arr, (math.log(nomi/dino))/(math.log(tAB[ip])))

print("\n q_arr looks like\n", q_arr)


#############################################################
###                  fiting the q_arr                   #####
#############################################################
def fit_q_f2(x, c0, c1, c2):
    delta_p = x - p_pcp_Greywall
    return c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f3(x, c0, c1, c2, c3):
    delta_p = x - p_pcp_Greywall
    return c3**(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0
 
def fit_q_f4(x, c0, c1, c2, c3, c4):
    delta_p = x - p_pcp_Greywall
    return c4*(delta_p**4) + c3**(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f5(x, c0, c1, c2, c3, c4, c5):
    delta_p = x - p_pcp_Greywall
    return c5*(delta_p**5) + c4*(delta_p**4) + c3**(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0

def fit_q_f6(x, c0, c1, c2, c3, c4, c5, c6):
    delta_p = x - p_pcp_Greywall
    return c6*(delta_p**6) + c5*(delta_p**5) + c4*(delta_p**4) + c3**(delta_p**3) + c2*(delta_p**2) + c1*(delta_p) + c0


popt2, pcov2 = cfit(fit_q_f2, pressure, q_arr)
popt3, pcov3 = cfit(fit_q_f3, pressure, q_arr)
popt4, pcov4 = cfit(fit_q_f4, pressure, q_arr)
popt5, pcov5 = cfit(fit_q_f5, pressure, q_arr)
popt6, pcov6 = cfit(fit_q_f6, pressure, q_arr)

print("\npopt2 looks like ",popt2,"\npopt3 looks like ",popt3,"\npopt4 looks like ",popt4,"\npopt5 looks like ",popt5,"\npopt6 looks like ",popt6)


#############################################################

fig, ax = plt.subplots(1,1)

#qfit2 = ax.plot(p_arr, fit_q_f2(p_arr, *popt2), 'b-.')
#qfit3 = ax.plot(p_arr, fit_q_f3(p_arr, *popt3), 'r-.')
qfit4 = ax.plot(p_arr, fit_q_f4(p_arr, *popt4), 'r:', label = r'$q_4 = {0:.6f} + {1:.6f}\Delta p + {2:.6f}\Delta p^2 + {3:.6f}\Delta p^3 + {4:.6}\Delta p^4$'.format(*popt4), linewidth = 4)
qfit5 = ax.plot(p_arr, fit_q_f5(p_arr, *popt5), 'm-', label = r'$q_5 = {0:.6f} + {1:.6f}\Delta p + {2:.6f}\Delta p^2 + {3:.6f}\Delta p^3 + {4:.6}\Delta p^4$ + {5:.6}\Delta p^5'.format(*popt5), linewidth = 4)
qfit6 = ax.plot(p_arr, fit_q_f6(p_arr, *popt6), 'c-.', label = r'$q_6 = {0:.6f} + {1:.6f}\Delta p + {2:.6f}\Delta p^2 + {3:.6f}\Delta p^3 + {4:.6}\Delta p^4 + {5:.6}\Delta p^5 + {5:.6f}\Delta p^6$'.format(*popt6), linewidth = 4)

# ax.plot(pressure, q_arr, linewidth = 2, color = "green")
q_expo = ax.scatter(pressure, q_arr, marker = "o", color = "green", s = 120, label = r"$q = ln(t_{AB}^{RWS29}/t_{AB}^{Greywall})$")
ax.set_xlabel(r"$p/bar$", fontsize = 16)
ax.set_ylabel(r"$q$", fontsize = 16)
ax.grid(True)
ax.legend(fontsize = 13.0)

# leg = ax.legend([qfit4],[r'$f4 = {0} + {1}\Delta p + {2}\Delta p^2 + {3}\Delta p^3 + {4}\Delta p^4$'.format(*popt4)], fontsize=15.0, loc='lower right')

#leg = ax.legend([qfit4, qfit5, qfit6, q_expo],[r'$f4 = {:.3g} + {:.3g}\Delta p + ({:.3g})\Delta p^2 + {:.3g}\Delta p^{3} + ({:.3g})\Delta p^4$'.format(*popt4)], fontsize=15.0, loc='lower right')


#    r"$$",r"$T_{c}^{PLTS}$",r"$p-T-IC, ConsP$",r"$p-T-HEC, ConsP$",r"$p-T-IC, ConsQ$",r"$p-T-HEC, ConsQ$"],  fontsize=15.0, loc='lower right')

plt.show()
