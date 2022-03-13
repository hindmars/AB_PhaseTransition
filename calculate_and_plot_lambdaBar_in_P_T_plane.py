############################################################
####        Important notations and descriptions       #####
############################################################


# This script is used for calculating and plotting the lambda_bar of the f(\phi) free energy
# deduced from original GL free energy by setting \phi \times D_{\alpha i} = A_{B} - A_{A}. Here A_{A} and A_{B}
# are particular order parameters of A and B phases. The explicit expression see the Mathematica note book.
# The resulted contour plot and density plot show \bar{\lambda} > 0.93 in whole p-T region with H = 0.

# comparing with Enqvist's 92' PRD, this result means Direct A-B tunneling free energy f(\phi) seems always belongs
# to " SSC " (Small Super Sooling) limit of Enqvist model.

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

###########################################################################################################
###########################################################################################################


# import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module
import Module_SC_Beta_V01 as SC_beta_gaps # all SC beta and gaps of A&B in SI unit
import Module_plot_TAB_line as TAB_line
import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math
from matplotlib import ticker, cm
import csv

# main code starts from here

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


###################################################
###          Temperature and pressure array    ####
###################################################

stepT = 0.02*(10**-3)*Kelvin 
Temperature = np.arange(0.0*(10**-3), 2.60*(10**-3), stepT) #Kelvin

stepPressure = 0.02*bar
pressure = np.arange(0.0*bar, 33.99*bar+stepPressure, stepPressure)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)

# saving array for lambdaBar
array_lambdaBar = np.zeros((lengthPressure,lengthT)) 

###################################################

for iP in range(0, lengthPressure, 1):
    print('\n\n Now P is:', pressure[iP], '\n\n')
    #indexT = math.floor(T/stepT)
    indexP = iP
    print('indexP is ',indexP)

    p = pressure[iP]

    for iT in range(0, lengthT, 1):
        
        indexT = iT
        print('indexT is',indexT, ' Temperature is ', Temperature[indexT])
    
        # call the SC_beta_gaps object
        alpha, beta1, beta2, beta3, beta4, beta5, betaA, betaB, DeltaA, DeltaB, phi0, Tcp, xitp, xi0p, t = SC_beta_gaps.calculate_SC_beta_gaps_etc(p,Temperature[indexT])

        # t = Temperature[indexT]/Tcp
        print('temperatureis:, ',t)
    
        if t >= 1:

           print(" bro, we just got temperature at Tc, save a np.nan. ")
      
           array_lambdaBar[indexP,indexT] = np.nan 

           
        else:

           # construct M2, delta, lambda 
           DeltaADeltaB_term = 8*(math.sqrt(2/3))*DeltaA*DeltaB
           betaAbetaB_coefficient = alpha/(3*betaA*betaB)

           # M2 
           M2 = (alpha/(2*phi0*phi0))*(DeltaADeltaB_term + betaAbetaB_coefficient*(14*(beta1 + beta2) + 7*beta3 + 3*beta4 + beta5))

           # delta
           delta_pT = ((3*alpha)/(2*phi0*phi0*phi0))*(DeltaADeltaB_term + betaAbetaB_coefficient*(8*beta1 + 14*beta2 + 5*beta3 + 7*beta4 + 5*beta5))

           # lambda
           lambda_pT = (alpha/(phi0*phi0*phi0*phi0))*(DeltaADeltaB_term + betaAbetaB_coefficient*(5*beta1 + 14*beta2 + 4*beta3 + 9*beta4 + 7*beta5))

           # lambda_bar
           lambda_bar = (9/2)*lambda_pT*(M2/(delta_pT*delta_pT))
           array_lambdaBar[indexP,indexT] = lambda_bar

           print(" \n lambda_bar is ", lambda_bar)

           

###########################################################################
################      contour plot lambda_bar   ###########################
###########################################################################

#LLLLL = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
LLLLL = [0.90, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.0, 1.005, 1.010, 1.015, 1.02, 1.025]

# LLLLL = 
X, Y = np.meshgrid(Temperature, pressure)
# fig, ax = plot1.subplots()
cs1 = plot1.contourf(X*(10**3), Y, array_lambdaBar, cmap=cm.PuBu_r, levels = LLLLL);
cs1a = plot1.contour(X*(10**3), Y, array_lambdaBar, levels = LLLLL, colors = 'yellow'); plot1.clabel(cs1a, inline=True, fontsize=9.5, colors='yellow')
plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
cb = plot1.colorbar(cs1)

# check the TAB line by using Module_TAB_line
cs2 = plot1.contour(X*(10**3), Y, array_lambdaBar, levels=[1.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=9.5, colors='k')
cs3 = plot1.contour(TAB_line.X*(10**3), TAB_line.Y, TAB_line.EnergyDensity_Difference_fABGL, levels=[0.0], colors='black'); plot1.clabel(cs3, inline=True, fontsize=9.5, colors='r')

# add title for the plot
plot1.title(r'$\bar{\lambda} = (9/2)({\lambda}M^{2}/{\delta}^{2})$ for $A \rightarrow B$. Red contour:$T_{AB}$')

plot1.savefig('lambdaBar_AB_Phase_Enqvist_model_1.pdf');
plot1.show()
           
plot1.clf()
plot1.cla()
plot1.close()

##############################################################################
###############         Density Plot of Lambda_bar     #######################
##############################################################################

DensityPlot = plot1.pcolormesh(X*(10**3), Y, array_lambdaBar);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
cb = plot1.colorbar(DensityPlot)

# plot the TAB line with lambda_bar = 1
cs2 = plot1.contour(X*(10**3), Y, array_lambdaBar, levels=[1.0], colors='red'); # plot1.clabel(cs2, inline=True, fontsize=9.5, colors='k')

# add title for the plot
plot1.title(r'$\bar{\lambda} = (9/2)({\lambda}M^{2}/{\delta}^{2})$ for $A \rightarrow B$, Red contour: $T_{AB}$')

plot1.savefig('lambdaBar_AB_Phase_Enqvist_model_1_DensityPlot.pdf');
plot1.show()


             
# plot1.scatter(n, normalized_ratio_fAfAB_for_LotynkExperimentData)
# plot1.scatter(n, normalized_ratio_ExperimentTension_for_LotynkExperimentData)
# plot1.scatter(n, normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)
# plot1.scatter(n, normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData)
# plot1.legend([r"$\frac{|f_{A}|}{f_{AB}}$" , r"$\frac{0.71 |f_{B}|}{f_{AB}}$", r"$E^{Thin-Wall}_{c}, \sigma_{AB} = 0.71|f_{B}|{\xi_{GL}}$", r"$E^{Thin-Wall}_{c}, \sigma_{AB} = |f_{A}|{\xi_{GL}}$"])
# plot1.show()


# ########################################################
# # contour plot energy barrier of the |f_A| xi evaluation
# ########################################################

# # LLLLL = [10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8, 3*10**8, 5*10**8, 7*10**8, 9*10**8, 10**9, 3*10**9, 5*10**9, 7*10**9, 9*10**9, 10**10]
# LLLLL = [10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8, 3*10**8, 5*10**8, 7*10**8, 9*10**8, 10**9, 3*10**9, 5*10**9, 7*10**9, 9*10**9, 10**10] 
# X, Y = np.meshgrid(Temperature, pressure)
# # fig, ax = plot1.subplots()
# cs1 = plot1.contourf(X*(10**3), Y, ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
# cb = plot1.colorbar(cs1)
# cb.set_ticks([10**5,5.0*10**5, 10**6, 5.0*10**6, 10**7, 5.0*10**7, 10**8, 5*10**8, 10**9, 5*10**9, 10**10])
# #cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi}$')
# # plot1.legend(loc='lower right')

# # Plot the TAB line
# cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

# # import the csv data of experment
# with open('slow_transition_lot_et_al_2021.csv', newline='') as f:
#     reader = csv.reader(f)
#     data = list(reader)

# print(data)


# list1 = list(zip(*data))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

# fuckData1 = []
# for ii in list(list1[0]):
#     fuckData1.append(float(ii))

# print("fuckData1 is", fuckData1)  
    
# fuckData2 = []
# for ii in list(list1[1]):
#     fuckData2.append(float(ii))    

# print("fuckData2 is", fuckData2)    

# plot1.scatter(fuckData1, fuckData2)
# #matplotlib.pyplot.scatter(fuckData1, fuckData2)

# # plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_{GL}})}^{3}}{{{\Delta}f_{AB}}^{2}} $')
# plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({|f_A|{\xi_{GL}}})^{3}}{{{\Delta}f_{AB}}^{2}} /k_{b}T $')
# plot1.savefig('Distribution_of_EnergyBarrier_ThinWall_Evaluation_Contour_fAxi_SCModule_versionI.pdf');

# plot1.show()

# ####################################################################
# # contour plot of energy barrier by 0.7 |f_B|xi tension
# ####################################################################


# X, Y = np.meshgrid(Temperature, pressure)
# # fig, ax = plot1.subplots()
# cs1 = plot1.contourf(X*(10**3), Y, ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
# cb = plot1.colorbar(cs1)
# cb.set_ticks([10**5,5.0*10**5, 10**6, 5.0*10**6, 10**7, 5.0*10**7, 10**8, 5*10**8, 10**9, 5*10**9, 10**10])
# #cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi}$')
# # plot1.legend(loc='lower right')

# cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

# # import the csv data of experment
# with open('slow_transition_lot_et_al_2021.csv', newline='') as f:
#     reader = csv.reader(f)
#     data = list(reader)

# print(data)


# list1 = list(zip(*data))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

# fuckData1 = []
# for ii in list(list1[0]):
#     fuckData1.append(float(ii))

# print("fuckData1 is", fuckData1)  
    
# fuckData2 = []
# for ii in list(list1[1]):
#     fuckData2.append(float(ii))    

# print("fuckData2 is", fuckData2)    

# plot1.scatter(fuckData1, fuckData2)
# #matplotlib.pyplot.scatter(fuckData1, fuckData2)

# plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_{GL}}})^{3}}{{{\Delta}f_{AB}}^{2}} /k_{b}T$')
# plot1.savefig('Distribution_of_EnergyBarrier_ThinWall_Evaluation_ExperimentalTension_Contour_SCModule_versionI.pdf');

# plot1.show()








# X, Y = np.meshgrid(Temperature, pressure)
# # cs1 = plot1.contourf(X*(10**3), Y, TempretureDependent_GL_CoherentLength*(10**9), locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
# DensityPlot = plot1.pcolormesh(X*(10**3), Y, TempretureDependent_GL_CoherentLength*(10**9));plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK')
# cb = plot1.colorbar(DensityPlot)
# plot1.clim(0, 250) # restrict the upper limit of colorbar
# cb.set_ticks([10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0]) # set the position of ticks of colorbar

# # contour plot
# LLLLL = np.arange(20, 200, 10)
# cs1 = plot1.contour(X*(10**3), Y, TempretureDependent_GL_CoherentLength*(10**9), locator=ticker.LogLocator(), levels=LLLLL, colors='black');plot1.clabel(cs1, inline=True, fontsize=9.5, colors='k')
# plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');

# # cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi}$')
# # plot1.legend(loc='lower right')

# # Plot the TAB line
# cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=9.5, colors='k')

# plot1.title(r'$\xi_{GL}(p, T)/nm$')
# plot1.savefig('Distribution_of_GL_CoherentLength_Contour_SCModule_versionI.pdf');

# plot1.show()



