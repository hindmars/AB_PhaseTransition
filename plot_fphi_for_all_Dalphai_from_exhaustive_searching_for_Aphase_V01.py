


import Module_SC_Beta_V01 as SC_beta_gaps # all SC beta and gaps of A&B in SI unit
# import Module_plot_TAB_line as TAB_line
# import Module_Lotynk_pT_parameter as Lotynk

import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math


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

fig, ax = plot1.subplots(1, 1)

for i in range(0, len(M2_filtered), 1):

    ax.plot((phi_array/DeltaA), fphi_saving_array[i, :])

plot1.ylim([-0.8, 2.5]);plot1.xlabel(r"$\phi/{\Delta_{A}}$");plot1.ylabel(r"$f(\phi)/J.m^{-3}$");ax.grid()

plot1.show()
    

    
####################################################################################################

# peak_filter = (phi_array/DeltaA) < 1.0
# print("\n\n peak_filter looks like : \n", peak_filter)
