


import Module_SC_Beta_V01 as SC_beta_gaps # all SC beta and gaps of A&B in SI unit
# import Module_plot_TAB_line as TAB_line
# import Module_Lotynk_pT_parameter as Lotynk

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from math import pi, exp
import math
import cmath


# unit in SI system

m = 1;s = 1; J = 1; Kelvin = 1; kg =1; bar = 1;# Length unit, Time unit, Energy unit, mass unit, pressure unit 
zeta3 = 1.2020569;
# kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1
kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1
c = 2.99792458*(10**(8)) # speed of light, m.s^-1

# hbar = 6.582119569*(10**(-16)) #plank constant, eV.s
hbar = 1.054571817*(10**(-34)) # planck constant, J.s
# u = 9.3149410242*(10**(8))*eV*(c**(-2)) # atomic mass unit, Dalton, eV.c^-2
u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom

#############################################################################################
#############   fucntion returning M2, \delta, lambda                   #####################
#############################################################################################

def calculator_of_lambdabar(d11, d12, d13, d21, d22, d23, d31, d32, d33,
                            psi11, psi12, psi13, psi21, psi22, psi23, psi31, psi32, psi33,
                            alpha, beta1, beta2, beta3, beta4, beta5,
                            gapA):
    
    M2_MonteCarlo = 2*(alpha - (-2*(beta2 + beta3 + beta4 + beta5) 
            + 2*beta3*(d13**2) 
            + 2*beta5*(d21**2 + d22**2 + d23**2 + d31**2 + d32**2 + d33**2)  
            + 2*beta1*(-1 + d13**2 + d21**2 + d22**2 + d23**2 + d31**2 + d32**2 + d33**2) 
            + (beta3 + beta4)*(d21**2 + d22**2 + 2*d23**2 + d31**2 + d32**2 + 2*d33**2))*(gapA**2) 
           + 2*(beta2 + beta4 + beta5)*(d11**2)*(gapA**2)*(np.cos(psi11))*(np.cos(psi11)) 
           + 4*(beta2 + beta4 + beta5)*d11*d12*(gapA**2)*(np.cos(psi11))*(np.sin(psi12))  
           - (gapA**2)*(-(beta5*(d21**2)*(np.cos(2*psi21))) 
                       + beta5*(d22**2)*(np.cos(2*psi22)) 
                       - beta5*(d31**2)*(np.cos(2*psi31)) 
                       + beta5*(d32**2)*(np.cos(2*psi32))
                       - 4*(beta1 + beta3)*d11*d12*(np.sin(psi11 - psi12)) 
                       - 2*(beta2 + beta4 + beta5)*(d12**2)*(np.sin(psi12))*(np.sin(psi12))  
                       + 2*d21*d22*((-beta3 + beta4)*(np.sin(psi21 - psi22)) - beta5*(np.sin(psi21 + psi22))) 
                       + 2*d31*d32*((-beta3 + beta4)*(np.sin(psi31 - psi32)) - beta5*(np.sin(psi31 + psi32)))))

    delta_MonteCarlo = 6*(np.sqrt(2))*gapA*(- (beta1*d11*(d22**2)*(np.cos(psi11 - 2*psi22)))  
                     - d12*d21*d22*(beta5*(np.cos(psi12 - psi21 - psi22)) + beta4*(np.cos(psi12 + psi21 - psi22)) + beta3*(np.cos(psi12 - psi21 + psi22)))  
                     - beta1*d11*(d23**2)*(np.cos(psi11 - 2*psi23)) 
                     - d13*d21*d23*(beta5*(np.cos(psi13 - psi21 - psi23)) + beta4*(np.cos(psi13 + psi21 - psi23)) + beta3*(np.cos(psi13 - psi21 + psi23))) 
                     - beta1*d11*(d32**2)*(np.cos(psi11 - 2*psi32)) 
                     - d12*d31*d32*(beta5*(np.cos(psi12 - psi31 - psi32)) + beta4*(np.cos(psi12 + psi31 - psi32)) + beta3*(np.cos(psi12 - psi31 + psi32)))  
                     - beta1*d11*(d33**2)*(np.cos(psi11 - 2*psi33)) 
                     - d13*d31*d33*(beta5*(np.cos(psi13 - psi31 - psi33)) + beta4*(np.cos(psi13 + psi31 - psi33)) + beta3*(np.cos(psi13 - psi31 + psi33))) 
                     + (np.cos(psi11))*(d11*(-beta2 + beta1*(-1 + d22**2 + d23**2 + d32**2 + d33**2)  
                                             +(beta3 + beta4 + beta5)*(-1 + d22**2 + d23**2 + d32**2 + d33**2)) 
                                        + 2*(beta1 + beta3)*(d12**3)*(np.sin(psi11 - psi12)))  
                     - (beta1 + beta3)*d12*(np.sin(2*psi11 - psi12)) 
                     - (beta2 + beta4 + beta5)*d12*(np.sin(psi12)) 
                     - 2*(beta1 + beta3)*d11*(d12**2)*(np.sin(psi11 - psi12))*(np.sin(psi12))  
                     + 2*(beta1 + beta3)*d12*(d13**2)*(np.cos(psi11 - psi12 + psi13))*(np.sin(psi11 - psi13)) 
                     - 2*(beta1 + beta3)*d11*(d13**2)*(np.sin(psi11 - psi13))*(np.sin(psi13))  
                     + d12*(d21**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (beta4 + beta5)*(np.sin(psi12)) + beta1*(np.sin(psi12 - 2*psi21)))  
                     - 2*(beta1 + beta5)*d11*(d21**2)*(np.sin(psi11 - psi21))*(np.sin(psi21))  
                     + d12*(d22**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (-beta3 + beta5)*(np.sin(psi12)) + (beta1 + beta5)*(np.sin(psi12 - 2*psi22)))  
                     + d11*d21*d22*(beta5*(np.sin(psi11 - psi21 - psi22)) - beta3*(np.sin(psi11 + psi21 - psi22)) - beta4*(np.sin(psi11 - psi21 + psi22)))  
                     + d12*(d23**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (beta4 + beta5)*(np.sin(psi12)) + beta1*(np.sin(psi12 - 2*psi23)))  
                     + d13*d22*d23*(beta5*(np.sin(psi13 - psi22 - psi23)) - beta4*(np.sin(psi13 + psi22 - psi23)) - beta3*(np.sin(psi13 - psi22 + psi23)))  
                     + d12*(d31**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (beta4 + beta5)*(np.sin(psi12)) + beta1*(np.sin(psi12 - 2*psi31)))  
                     - 2*(beta1 + beta5)*d11*(d31**2)*(np.sin(psi11 - psi31))*(np.sin(psi31))  
                     + d12*(d32**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (-beta3 + beta5)*(np.sin(psi12)) + (beta1 + beta5)*(np.sin(psi12 - 2*psi32)))  
                     + d11*d31*d32*(beta5*(np.sin(psi11 - psi31 - psi32)) - beta3*(np.sin(psi11 + psi31 - psi32)) - beta4*(np.sin(psi11 - psi31 + psi32)))  
                     + d12*(d33**2)*((beta1 + beta3)*(np.sin(2*psi11 - psi12)) + (beta4 + beta5)*(np.sin(psi12)) + beta1*(np.sin(psi12 - 2*psi33)))  
                     + d13*d32*d33*(beta5*(np.sin(psi13 - psi32 - psi33)) - beta4*(np.sin(psi13 + psi32 - psi33)) - beta3*(np.sin(psi13 - psi32 + psi33))))

    delta2_MonteCarlo = np.square(delta_MonteCarlo)

    lambda_MonteCarlo = 4*(beta2 + beta4 
                        + (beta1 + beta3 + beta5)*(d11**4) + (beta1 + beta3 + beta5)*(d12**4) + (beta1 + beta3 + beta5)*(d13**4) 
                        - 2*beta4*(d13**2)*(d21**2) + (beta1 + beta3 + beta5)*(d21**4) 
                        - 2*beta4*(d22**2) + 2*(beta4 + beta5)*(d21**2)*(d22**2) 
                        + (beta1 + beta3 + 2*beta4 + beta5)*(d22**4) 
                        - 2*beta4*(d23**2) + 2*(beta3 + beta4)*(d13**2)*(d23**2) 
                        + 2*(beta4 + beta5)*(d21**2)*(d23**2) + 2*(2*beta4 + beta5)*(d22**2)*(d23**2) 
                        + (beta1 + beta3 + 2*beta4 + beta5)*(d23**4) 
                        - 2*beta4*(d13**2)*(d31**2) + 2*beta3*(d21**2)*(d31**2) + (beta1 + beta3 + beta5)*(d31**4) + 2*(d11**2)*(beta5*(d12**2 + d13**2) + beta3*(d21**2 + d31**2)) 
                        - 2*beta4*(d32**2) + 2*(beta3 + 2*beta4)*(d22**2)*(d32**2) + 2*beta4*(d23**2)*(d32**2) + 2*(beta4 + beta5)*(d31**2)*(d32**2) + (beta1 + beta3 + 2*beta4 + beta5)*(d32**4) 
                        + 2*(d12**2)*(beta5*(d13**2) + (beta3 + beta4)*(d22**2) - beta4*(d21**2 + d31**2) + (beta3 + beta4)*(d32**2)) 
                        - 2*beta4*(d33**2) + 2*(beta3 + beta4)*(d13**2)*(d33**2) + 2*beta4*(d22**2)*(d33**2) + 2*(beta3 + 2*beta4)*(d23**2)*(d33**2) 
                        + 2*(beta4 + beta5)*(d31**2)*(d33**2) + 2*(2*beta4 + beta5)*(d32**2)*(d33**2) 
                        + (beta1 + beta3 + 2*beta4 + beta5)*(d33**4) 
                        + 2*((beta1 + beta3)*(d11**2)*(d12**2)*(np.cos(2*psi11 - 2*psi12)) + (beta1 + beta3)*(d11**2)*(d13**2)*(np.cos(2*psi11 - 2*psi13)) 
                             + beta1*(d12**2)*(d13**2)*(np.cos(2*psi12 - 2*psi13)) + beta3*(d12**2)*(d13**2)*(np.cos(2*psi12 - 2*psi13)) + beta1*(d11**2)*(d21**2)*(np.cos(2*psi11 - 2*psi21)) 
                             + beta5*(d11**2)*(d21**2)*(np.cos(2*psi11 - 2*psi21)) + beta1*(d12**2)*(d21**2)*(np.cos(2*psi12 - 2*psi21)) + beta1*(d13**2)*(d21**2)*(np.cos(2*psi13 - 2*psi21)) 
                             + beta1*(d11**2)*(d22**2)*(np.cos(2*psi11 - 2*psi22)) + beta1*(d12**2)*(d22**2)*(np.cos(2*psi12 - 2*psi22)) + beta5*(d12**2)*(d22**2)*(np.cos(2*psi12 - 2*psi22)) 
                             + beta1*(d13**2)*(d22**2)*(np.cos(2*psi13 - 2*psi22)) + beta1*(d21**2)*(d22**2)*(np.cos(2*psi21 - 2*psi22)) + beta3*(d21**2)*(d22**2)*(np.cos(2*psi21 - 2*psi22)) 
                             + 2*beta5*d11*d12*d21*d22*(np.cos(psi11 + psi12 - psi21 - psi22)) + 2*beta3*d11*d12*d21*d22*(np.cos(psi11 - psi12 + psi21 - psi22)) 
                             + 2*beta4*d11*d12*d21*d22*(np.cos(psi11 - psi12 - psi21 + psi22)) + beta1*(d11**2)*(d23**2)*(np.cos(2*psi11 - 2*psi23)) + beta1*(d12**2)*(d23**2)*(np.cos(2*psi12 - 2*psi23)) 
                             + beta1*(d13**2)*(d23**2)*(np.cos(2*psi13 - 2*psi23)) + beta5*(d13**2)*(d23**2)*(np.cos(2*psi13 - 2*psi23)) + beta1*(d21**2)*(d23**2)*(np.cos(2*psi21 - 2*psi23)) 
                             + beta3*(d21**2)*(d23**2)*(np.cos(2*psi21 - 2*psi23)) + beta1*(d22**2)*(d23**2)*(np.cos(2*psi22 - 2*psi23)) + beta3*(d22**2)*(d23**2)*(np.cos(2*psi22 - 2*psi23)) 
                             + 2*beta5*d11*d13*d21*d23*(np.cos(psi11 + psi13 - psi21 - psi23)) + 2*beta3*d11*d13*d21*d23*(np.cos(psi11 - psi13 + psi21 - psi23)) 
                             + 2*beta5*d12*d13*d22*d23*(np.cos(psi12 + psi13 - psi22 - psi23)) + 2*beta3*d12*d13*d22*d23*(np.cos(psi12 - psi13 + psi22 - psi23)) 
                             + 2*beta4*d11*d13*d21*d23*(np.cos(psi11 - psi13 - psi21 + psi23)) + 2*beta4*d12*d13*d22*d23*(np.cos(psi12 - psi13 - psi22 + psi23)) 
                             + beta1*(d11**2)*(d31**2)*(np.cos(2*psi11 - 2*psi31)) + beta5*(d11**2)*(d31**2)*(np.cos(2*psi11 - 2*psi31)) + beta1*(d12**2)*(d31**2)*(np.cos(2*psi12 - 2*psi31)) 
                             + beta1*(d13**2)*(d31**2)*(np.cos(2*psi13 - 2*psi31)) + beta1*(d21**2)*(d31**2)*(np.cos(2*psi21 - 2*psi31)) + beta5*(d21**2)*(d31**2)*(np.cos(2*psi21 - 2*psi31)) 
                             + beta1*(d22**2)*(d31**2)*(np.cos(2*psi22 - 2*psi31)) + beta1*(d23**2)*(d31**2)*(np.cos(2*psi23 - 2*psi31)) + beta1*(d11**2)*(d32**2)*(np.cos(2*psi11 - 2*psi32)) 
                             + beta1*(d12**2)*(d32**2)*(np.cos(2*psi12 - 2*psi32)) + beta5*(d12**2)*(d32**2)*(np.cos(2*psi12 - 2*psi32)) + beta1*(d13**2)*(d32**2)*(np.cos(2*psi13 - 2*psi32)) 
                             + beta1*(d21**2)*(d32**2)*(np.cos(2*psi21 - 2*psi32)) + beta1*(d22**2)*(d32**2)*(np.cos(2*psi22 - 2*psi32)) + beta5*(d22**2)*(d32**2)*(np.cos(2*psi22 - 2*psi32)) 
                             + beta1*(d23**2)*(d32**2)*(np.cos(2*psi23 - 2*psi32)) + beta1*(d31**2)*(d32**2)*(np.cos(2*psi31 - 2*psi32)) + beta3*(d31**2)*(d32**2)*(np.cos(2*psi31 - 2*psi32)) 
                             + 2*beta5*d11*d12*d31*d32*(np.cos(psi11 + psi12 - psi31 - psi32)) + 2*beta5*d21*d22*d31*d32*(np.cos(psi21 + psi22 - psi31 - psi32)) 
                             + 2*beta3*d11*d12*d31*d32*(np.cos(psi11 - psi12 + psi31 - psi32)) + 2*beta3*d21*d22*d31*d32*(np.cos(psi21 - psi22 + psi31 - psi32)) 
                             + 2*beta4*d11*d12*d31*d32*(np.cos(psi11 - psi12 - psi31 + psi32)) + 2*beta4*d21*d22*d31*d32*(np.cos(psi21 - psi22 - psi31 + psi32)) 
                             + beta1*(d11**2)*(d33**2)*(np.cos(2*psi11 - 2*psi33)) + beta1*(d12**2)*(d33**2)*(np.cos(2*psi12 - 2*psi33)) + beta1*(d13**2)*(d33**2)*(np.cos(2*psi13 - 2*psi33)) 
                             + beta5*(d13**2)*(d33**2)*(np.cos(2*psi13 - 2*psi33)) + beta1*(d21**2)*(d33**2)*(np.cos(2*psi21 - 2*psi33)) + beta1*(d22**2)*(d33**2)*(np.cos(2*psi22 - 2*psi33)) 
                             + beta1*(d23**2)*(d33**2)*(np.cos(2*psi23 - 2*psi33)) + beta5*(d23**2)*(d33**2)*(np.cos(2*psi23 - 2*psi33)) + beta1*(d31**2)*(d33**2)*(np.cos(2*psi31 - 2*psi33)) 
                             + beta3*(d31**2)*(d33**2)*(np.cos(2*psi31 - 2*psi33)) + beta1*(d32**2)*(d33**2)*(np.cos(2*psi32 - 2*psi33)) + beta3*(d32**2)*(d33**2)*(np.cos(2*psi32 - 2*psi33)) 
                             + 2*beta5*d11*d13*d31*d33*(np.cos(psi11 + psi13 - psi31 - psi33)) + 2*beta5*d21*d23*d31*d33*(np.cos(psi21 + psi23 - psi31 - psi33)) 
                             + 2*beta3*d11*d13*d31*d33*(np.cos(psi11 - psi13 + psi31 - psi33)) + 2*beta3*d21*d23*d31*d33*(np.cos(psi21 - psi23 + psi31 - psi33)) 
                             + 2*beta5*d12*d13*d32*d33*(np.cos(psi12 + psi13 - psi32 - psi33)) + 2*beta5*d22*d23*d32*d33*(np.cos(psi22 + psi23 - psi32 - psi33)) 
                             + 2*beta3*d12*d13*d32*d33*(np.cos(psi12 - psi13 + psi32 - psi33)) + 2*beta3*d22*d23*d32*d33*(np.cos(psi22 - psi23 + psi32 - psi33)) 
                             + 2*beta4*d11*d13*d31*d33*(np.cos(psi11 - psi13 - psi31 + psi33)) + 2*beta4*d21*d23*d31*d33*(np.cos(psi21 - psi23 - psi31 + psi33)) 
                             + 2*beta4*d12*d13*d32*d33*(np.cos(psi12 - psi13 - psi32 + psi33)) + 2*beta4*d22*d23*d32*d33*(np.cos(psi22 - psi23 - psi32 + psi33))))

    lambda_bar_MonteCarlo = (9/2)*lambda_MonteCarlo*(M2_MonteCarlo / delta2_MonteCarlo)

    return M2_MonteCarlo, delta_MonteCarlo, lambda_MonteCarlo, lambda_bar_MonteCarlo

##############################################################################################
######## load the elements of D_alphai getting from exhaustive searching of A-phase ##########
##############################################################################################

with np.load('filtered_lambdaBar_Dalphai_20000loops.npz') as data:
    d11 = data['d11']
    d12 = data['d12']
    d13 = data['d13']
    d21 = data['d21']
    d22 = data['d22']
    d23 = data['d23']
    d31 = data['d31']
    d32 = data['d32']
    d33 = data['d33']

    psi11 = data['psi11']
    psi12 = data['psi12']
    psi13 = data['psi13']
    psi21 = data['psi21']
    psi22 = data['psi22']
    psi23 = data['psi23']
    psi31 = data['psi31']
    psi32 = data['psi32']
    psi33 = data['psi33']
    
    lambda_bar = data['lambdaBar']

    M2 = data['M2']
    delta = data['delta']
    lambda_pT = data['lambda_pT']
    
    
print(" the number of nagative elements of \delta is : ", (delta[delta < 0]).size)

###############################################################################################

###############################################################################################
##################             positive lambda filter            ##############################
###############################################################################################

positive_lambda_filter = delta > 0
print(" \n\n the positive lambda filter is : \n ", positive_lambda_filter)

lambda_bar_filtered = lambda_bar[positive_lambda_filter]

d11_filtered = d11[positive_lambda_filter]
d12_filtered = d12[positive_lambda_filter]
d13_filtered = d13[positive_lambda_filter]
d21_filtered = d21[positive_lambda_filter]
d22_filtered = d22[positive_lambda_filter]
d23_filtered = d23[positive_lambda_filter]
d31_filtered = d31[positive_lambda_filter]
d32_filtered = d32[positive_lambda_filter]
d33_filtered = d33[positive_lambda_filter]


psi11_filtered = psi11[positive_lambda_filter]
psi12_filtered = psi12[positive_lambda_filter]
psi13_filtered = psi13[positive_lambda_filter]
psi21_filtered = psi21[positive_lambda_filter]
psi22_filtered = psi22[positive_lambda_filter]
psi23_filtered = psi23[positive_lambda_filter]
psi31_filtered = psi31[positive_lambda_filter]
psi32_filtered = psi32[positive_lambda_filter]
psi33_filtered = psi33[positive_lambda_filter]


M2_filtered = M2[positive_lambda_filter]
delta_filtered = delta[positive_lambda_filter]
lambda_pT_filtered = lambda_pT[positive_lambda_filter]


##################################################################################################

##################################################################################################
####################      prepare the \phi for p = 32bar, T = 0.1 mK    ##########################
##################################################################################################

p = 32.0*bar
Temperature = 0.1*(10**(-3))*Kelvin

alpha, beta1, beta2, beta3, beta4, beta5, betaA, betaB, DeltaA, DeltaB, phi0, Tcp, xitp, xi0p, t = SC_beta_gaps.calculate_SC_beta_gaps_etc(p,Temperature)

print(" \n\n gapA is : ", DeltaA)

# np.array of \phi
phi_array = np.arange(0.0,2.0,0.01) * DeltaA

print(" \n\n phi_array looks like : \n", phi_array)

##################################################################################################

##################################################################################################
##########   calculate out all fphi for all D_alphaii, with M2, delta and lambda_pT  #############
##################################################################################################

print(" \n\n len(lambda_bar_filtered) is :", len(lambda_bar_filtered), " len(phi_array) is :", len(phi_array))
fphi_saving_array = np.zeros((0,len(phi_array)), dtype = np.float64)
# fphi_saving_array = np.array([])

for i in range(0, len(M2_filtered), 1):

    print(" \n now i is : ", i)


    fphi = (1.0/2.0) * M2_filtered[i] * (phi_array**2) - (1.0/3.0) * delta_filtered[i] * (phi_array**3) + (1.0/4.0) * lambda_pT_filtered[i] * (phi_array**4)

    # print(" \n\n fphi array looks like : ", fphi)

    fphi_saving_array = np.append(fphi_saving_array, [fphi], axis=0)


# print(" \n\n shape of fphi_saving_array : ", fphi_saving_array.shape)

# print(" \n\n piece of fphi_saving_array : ", fphi_saving_array[0:5, :])




###################################################################################################

###################################################################################################
############         plot all fphi in every row of fphi_saving_array                ###############
###################################################################################################

# print(" the oth row of fphi_saving_array looks like : \ ", fphi_saving_array[0, :], " \n the 100th row of fphi_saving_array looks like : \n ", fphi_saving_array[100, :])

# fig, ax = plt.subplots(1, 1)

# for i in range(0, len(M2_filtered), 1):

#     ax.plot((phi_array/DeltaA), fphi_saving_array[i, :])

# plt.ylim([-0.8, 2.5]);plt.xlabel(r"$\phi/{\Delta_{A}}$");plt.ylabel(r"$f(\phi)/J.m^{-3}$");ax.grid()

# plt.show()
    

    
####################################################################################################

# peak_filter = (phi_array/DeltaA) < 1.0
# print("\n\n peak_filter looks like : \n", peak_filter)

# cutted_max = np.array([])

# for i in range(0, len(M2_filtered), 1):

#     fphi_i = fphi_saving_array[i, :]

#     cutted_fphi = fphi_i[peak_filter]

#     print("\n max of this cutted fphi is :\n", cutted_fphi.max(axis=0))

#     cutted_max = np.append(cutted_max, cutted_fphi.max(axis=0))

#     print("\n\n now fphi_i is :\n", fphi_i, " \n now cutted_fphi is : \n ", cutted_fphi)


# fig1, ax1 = plt.subplots(1,1)    
# ax1.scatter(np.arange(0, len(cutted_max), 1), cutted_max);ax1.grid()
# # plt.show()


# fig, ax = plt.subplots(1, 1)

# ax.plt((phi_array/DeltaA), fphi_saving_array[1307, :], (phi_array/DeltaA), fphi_saving_array[3929, :]);ax.grid()
# plt.ylim([-0.8, 0.5]);
# plt.xlim([0, 1.0]);
# plt.show()


# let's show out the D_alphai with lowest barrier

# index = 1307
# Dalphaii1307 = np.array([ [d11_filtered[index]*cmath.exp(psi11_filtered[index]*1j), d12_filtered[index]*cmath.exp(psi12_filtered[index]*1j), d13_filtered[index]*cmath.exp(psi13_filtered[index]*1j)],
#                           [d21_filtered[index]*cmath.exp(psi21_filtered[index]*1j), d22_filtered[index]*cmath.exp(psi22_filtered[index]*1j), d23_filtered[index]*cmath.exp(psi23_filtered[index]*1j)],
#                           [d31_filtered[index]*cmath.exp(psi31_filtered[index]*1j), d32_filtered[index]*cmath.exp(psi32_filtered[index]*1j), d33_filtered[index]*cmath.exp(psi33_filtered[index]*1j)]
#                         ])
# print(" \n\n the 1307th Dalphai looks like : \n ", Dalphaii1307)

# index = 3929
# Dalphaii3929 = np.array([ [d11_filtered[index]*cmath.exp(psi11_filtered[index]*1j), d12_filtered[index]*cmath.exp(psi12_filtered[index]*1j), d13_filtered[index]*cmath.exp(psi13_filtered[index]*1j)],
#                           [d21_filtered[index]*cmath.exp(psi21_filtered[index]*1j), d22_filtered[index]*cmath.exp(psi22_filtered[index]*1j), d23_filtered[index]*cmath.exp(psi23_filtered[index]*1j)],
#                           [d31_filtered[index]*cmath.exp(psi31_filtered[index]*1j), d32_filtered[index]*cmath.exp(psi32_filtered[index]*1j), d33_filtered[index]*cmath.exp(psi33_filtered[index]*1j)]
#                         ])
# print(" \n\n the 3929th Dalphai looks like : \n ", Dalphaii3929)

p = 32.0*bar
Temperature = 0.1*(10**(-3))*Kelvin

index = 3929

alpha, beta1, beta2, beta3, beta4, beta5, betaA, betaB, DeltaA, DeltaB, phi0, Tcp, xitp, xi0p, t = SC_beta_gaps.calculate_SC_beta_gaps_etc(p,Temperature)

M2_3929, delta_3929, lambda_3929, lambda_bar_3929 = calculator_of_lambdabar(d11_filtered[index], d12_filtered[index], d13_filtered[index],
                                                                                d21_filtered[index], d22_filtered[index], d23_filtered[index],
                                                                                d31_filtered[index], d32_filtered[index], d33_filtered[index],
                                                                                psi11_filtered[index], psi12_filtered[index], psi13_filtered[index],
                                                                                psi21_filtered[index], psi22_filtered[index], psi23_filtered[index],
                                                                                psi31_filtered[index], psi32_filtered[index], psi33_filtered[index],
                                                                                alpha, beta1, beta2, beta3, beta4, beta5,
                                                                                DeltaA)

print(" \n\n the 3929th D_alphai gives out : M2 =", M2_3929, " \delta = ", delta_3929," lambda = ", lambda_3929, " lambda_bar = ",lambda_bar_3929)

print(" \n\n the No. 3929 elements in M2_filtered is ", M2_filtered[index], " \delta is ", delta_filtered[index], " \lambda is ", lambda_pT_filtered[index], " lambda_bar is ", lambda_bar_filtered[index])

fig, ax = plt.subplots(1, 1)
fphi3929 = (1.0/2.0) * M2_3929 * (phi_array**2) - (1.0/3.0) * delta_3929 * (phi_array**3) + (1.0/4.0) * lambda_3929 * (phi_array**4)

ax.plot(phi_array/DeltaA, fphi3929);
ax.set_xlabel(r"$\phi/{\Delta_{A}}$");ax.set_ylabel(r"$f(\phi)/J.m^{-3}$");ax.grid(True)
ax.set_ylim([-0.5, 1.0]);ax.set_xlim([0, 1.2]);ax.legend(labels=['lowest barrier line'], loc = 'upper left')

plt.show()
