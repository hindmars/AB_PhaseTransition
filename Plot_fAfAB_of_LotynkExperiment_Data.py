#######################
# Important notations #
#######################


# This script is used for calculating and plotting the thin-wall energy barriers and |f_A|/f_AB
# under Lotynk's experimental pressure and temperature. Mark guessed the " accidence " between
# |f_A|/f_AB plot and experimental nucleation events is clue about some kinds of " thick wall " bubble
# in Enqvist's 1992 paper.

# Then we actually expected the Thin-Wall energy barriers have no relationship with |f_A|/f_AB, which is monotonically
# determine the 3D action S_3 of " thick wall " bubble in Enqvist's paper. And the resulted scatter plot exactly shows
# out thin-wall energy barrier has poor relation with |f_A|/f_AB. and the former has bigger stand deviation (0.3037) than
# the later (0.1). If the probalility of neucleations are same, then this bigger deviation just means the bubble definatley
# is not thin-wall bubble.

##############################################################################################################

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

# This script uses the SC module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version.2 of ebergy difference code, in which the pico-Joule(pJ 10^-12) unit is used 

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569;


import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module
import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math
from matplotlib import ticker, cm
import csv

# main code starts from here

m = 1;s = 1; J = 1; Kelvin = 1 # Length unit, Time unit, Energy unit
bar = 1
zeta3 = 1.2020569;
# kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1
kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1
c = 2.99792458*(10**(8)) # speed of light, m.s^-1

# hbar = 6.582119569*(10**(-16)) #plank constant, eV.s
hbar = 1.054571817*(10**(-34)) # planck constant, J.s
# u = 9.3149410242*(10**(8))*eV*(c**(-2)) # atomic mass unit, Dalton, eV.c^-2
u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom


###########################################
# build object, and check the intepolatinon 

BetaObject = SC.BETA('betaAndTc')

#p = 21.1 # pressure, bar
# pressureArray = np.arange(0.0, 34.0, 2.0)

#for p in pressureArray:
    
###################################################

# import the csv data of experment
with open('slow_transition_lot_et_al_2021.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

print(data)


list1 = list(zip(*data))
print(list1)
print(list(list1[0]))
print(list1[1])

fuckData1 = []
for ii in list(list1[0]):
    fuckData1.append(float(ii))

print("fuckData1 is", fuckData1)  
    
fuckData2 = []
for ii in list(list1[1]):
    fuckData2.append(float(ii))    

print("fuckData2 is", fuckData2)    

# plot1.scatter(fuckData1, fuckData2)

# plot1.show()

# calculate free energy under different temperature


# stepT = 0.001*(10**-3) 
# Temperature = np.arange(2.0*(10**-3), 2.40*(10**-3), stepT) #Kelvin
Temperature = np.asarray(fuckData1)*(10**(-3)) # Kelvin

# stepPressure = 0.001*bar
# # pressure = np.arange(20.0, 33.9*bar, stepPressure)
# pressure = np.arange(19.9, 22.99*bar+stepPressure, stepPressure)
pressure = np.asarray(fuckData2)*bar

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)




###################################################
# saving arraies for Lotynk's experiment parameters
###################################################

ratio_fAfAB_for_LotynkExperimentData = np.zeros((1,lengthPressure))

ratio_ExperimentTension_for_LotynkExperimentData = np.zeros((1,lengthPressure))

ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData = np.zeros((1,lengthPressure))

ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData = np.zeros((1,lengthPressure))

###################################################

for iP in range(0, lengthPressure, 1):
    print('\n\n Now P is:', pressure[iP], '\n\n')
    #indexT = math.floor(T/stepT)
    indexP = iP
    print('indexP is ',indexP)

    p = pressure[iP]
    
    BetaObject.c1_function(SC.P,SC.c1,p);c1p = BetaObject.c1p
    BetaObject.c2_function(SC.P,SC.c2,p);c2p = BetaObject.c2p
    BetaObject.c3_function(SC.P,SC.c3,p);c3p = BetaObject.c3p
    BetaObject.c4_function(SC.P,SC.c4,p);c4p = BetaObject.c4p
    BetaObject.c5_function(SC.P,SC.c5,p);c5p = BetaObject.c5p
    BetaObject.tc_function(SC.P,SC.Tc,p);Tcp = (BetaObject.tcp)*(10**(-3)) # turn mK to Kelvin 10^(-3)
    BetaObject.mStar_function(SC.P,SC.Ms,p); mEffective = (BetaObject.ms)*m3
    BetaObject.vFermi_function(SC.P,SC.VF,p); vFermi = BetaObject.vf
    BetaObject.xi0_function(SC.P,SC.XI0,p); xi0p = (BetaObject.xi0)*(10**(-9)) # turn nm to m 10^(-9)
    

    c245p = c2p + c4p + c5p;c12p = c1p + c2p;c345p = c3p + c4p + c5p

    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,' ,xi0p is ', xi0p, '\n\n')

    N0 = ((mEffective**(2))*vFermi)/((2*pi*pi)*(hbar**(3))) # energy density of Fermi surface

    print('\npressure is, ',p,' effective mass is, ', mEffective, ' Fermi velocity is,', vFermi, ' N(0) is ',N0,'\n\n')
    
#    for iT in range(0, lengthT, 1):
        #indexDelta = math.floor(delta/stepDelta)
#        indexT = iT
#        print('indexT is',indexT, ' Temperature is ', Temperature[indexT])

    t = Temperature[iP]/Tcp
    print('temperatureis:, ',t)

    xiGLWeakCoupling = math.sqrt((7.0*zeta3)/20.0)*xi0p # weak coupling GL coherent length in PRB 101 024517

         

    if t >= 1:

      print(" bro, we just got temperature at Tc, save a np.nan. ")
      
      ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData[0,indexP] = np.nan

      ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData[0,indexP] = np.nan

      ratio_fAfAB_for_LotynkExperimentData[0,indexP] = np.nan

      ratio_ExperimentTension_for_LotynkExperimentData[0,indexP] = np.nan
      

           
    else:

        xitp = xiGLWeakCoupling/math.sqrt(1-t) # temperature dependent coherent length  

        
           
        print(' Temperature dependent GL xi is ', xitp)
            
        alphaRed = (1/3)*(t-1)

        beta245Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(2+t*c245p) # A Phase
        beta12345Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(3.0*(1+t*c12p) + (2+t*c345p)) # B phase

        deltaA = math.sqrt((-alphaRed)/(4*beta245Red)) # A -> B

        deltaB = math.sqrt((-alphaRed)/(2*beta12345Red)) # A -> B
        
        fAGLRed = alphaRed*2*(deltaA**2) + 4*beta245Red*(deltaA**4)
 
#       fBGLRed = alphaRed*3*(delta**2) + 3*beta12345Red*(delta**4)
        fBGLRed = alphaRed*3*(deltaB**2) + 3*beta12345Red*(deltaB**4)
           
        fAGL = fAGLRed * N0 # real value of fAGL
        fBGL = fBGLRed * N0 # real value of fBGL
        DiffFABGL = fAGL - fBGL # real value of \Delta f_AB

             
           
        print('fAGL in SI unit is:', fAGL)
        print('fBGL in SI unit is:', fBGL)
        print('DiffFABGL in SI unit is:', DiffFABGL)

           
           #######################################################################
           # calculation for thin-wall evaluation of energy barrier in unit of KbT
           #######################################################################
        numerator1 = (abs(fAGL)*xitp)**3 # tension for |fA| xi evaluation
        numerator2 = (0.71*abs(fBGL)*xitp)**3 # tension from experiment closed TAB
        dinorminator = DiffFABGL**2 
        KbT = kb*Temperature[iP]
        print('KbT is:', KbT)

        EnergyBarrier_fAxi = (16/3)*pi*(numerator1/dinorminator) # Energy Barrier with SI unit, formule see 25th Nov. 2021 note
        EnergyBarrier_experimentTension = (16/3)*pi*(numerator2/dinorminator) 

        print(' energy barrier by |fA|xi in SI unit is ', EnergyBarrier_fAxi)
        print(' energy barrier by 0.7 |f_B|xi in SI unit is ', EnergyBarrier_experimentTension)
           
        EnergyBarrier_fAxi_KbT = EnergyBarrier_fAxi/KbT # Energy Barrier with SI unit
        print(' energy barrier by |fA|xi in unit of KbT is ', EnergyBarrier_fAxi_KbT)

        EnergyBarrier_experimentTension_KbT = EnergyBarrier_experimentTension/KbT # Energy Barrier with SI unit
        print(' energy barrier by 0.7 |fB| xi in unit of KbT is ', EnergyBarrier_experimentTension_KbT)

           
           ###################################################################
           # calculation for |f_A|/f_AB and 0.7 |f_B|/f_AB
           ###################################################################
           
        ratio_fAfAB = abs(fAGL)/DiffFABGL

        ratio_ExperimentTension = (0.7*abs(fBGL))/DiffFABGL # 0.7 |f_B| \xiGL /(f_AB \xi_GL)

           ###################################################################
        
        if DiffFABGL > 0:

              
            ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData[0,indexP] = EnergyBarrier_fAxi_KbT # in unit of |f_B| \xi_GL^3

            print(" thin wall energy barrier in KbT is ", ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData)

            ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData[0,indexP] = EnergyBarrier_experimentTension_KbT # in unit of |f_B| \xi_GL^3

            print(" thin wall energy barrier in KbT is ", ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)

              
              

            ratio_fAfAB_for_LotynkExperimentData[0,indexP] = ratio_fAfAB
              
            print(" ratio |f_A|/f_AB is ", ratio_fAfAB_for_LotynkExperimentData)

            ratio_ExperimentTension_for_LotynkExperimentData[0,indexP] = ratio_ExperimentTension

            print(" ratio 0.7 |f_B|/f_AB is ", ratio_ExperimentTension_for_LotynkExperimentData)


        else:

            print(" bro, you get negative logarithm, save a np.nan. ")

            ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData[0,indexP] = np.nan

            ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData[0,indexP] = np.nan

            ratio_fAfAB_for_LotynkExperimentData[0,indexP] = np.nan

            ratio_ExperimentTension_for_LotynkExperimentData[0,indexP] = np.nan
        
print(" \n\n\n the mean of |f_A|/f_AB is : ", ratio_fAfAB_for_LotynkExperimentData.mean())
print("  the variance of |f_A|/f_AB is : ", ratio_fAfAB_for_LotynkExperimentData.var())
print("  the stand diviation of |f_A|/f_AB is : ", ratio_fAfAB_for_LotynkExperimentData.std())

normalized_ratio_fAfAB_for_LotynkExperimentData = ratio_fAfAB_for_LotynkExperimentData/(ratio_fAfAB_for_LotynkExperimentData.mean())
print(" \n\n\n the normalized |f_A|/f_AB is : ", normalized_ratio_fAfAB_for_LotynkExperimentData)
print(" \nthe variance of normalized |f_A|/f_AB is : ", normalized_ratio_fAfAB_for_LotynkExperimentData.var())
print(" the stand diviation of normalized |f_A|/f_AB is : ", normalized_ratio_fAfAB_for_LotynkExperimentData.std())

##################################################################################


print(" \n\n\n the mean of 0.71 |f_B|/f_AB is : ", ratio_ExperimentTension_for_LotynkExperimentData.mean())
print("  the variance of 0.71 |f_B|/f_AB is : ", ratio_ExperimentTension_for_LotynkExperimentData.var())
print("  the stand diviation of 0.71 |f_B|/f_AB is : ", ratio_ExperimentTension_for_LotynkExperimentData.std())

normalized_ratio_ExperimentTension_for_LotynkExperimentData = ratio_ExperimentTension_for_LotynkExperimentData/(ratio_ExperimentTension_for_LotynkExperimentData.mean())
print(" \n\n\n the normalized 0.71 |f_B|/f_AB is : ", normalized_ratio_ExperimentTension_for_LotynkExperimentData)
print(" \nthe variance of normalized 0.71 |f_B|/f_AB is : ", normalized_ratio_ExperimentTension_for_LotynkExperimentData.var())
print(" the stand diviation of normalized 0.71 |f_B|/f_AB is : ", normalized_ratio_ExperimentTension_for_LotynkExperimentData.std())

##################################################################################

print(" \n\n\n the mean of Thin-Wall Energy Barrier |f_A|/f_ABxi is : ", ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.mean())
print("  the variance of Thin-Wall Energy Barrier(experiemt) is : ", ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.var())
print("  the stand diviation of Thin-Wall Energy Barrier(experiemt) is : ", ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.std())

normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData = ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData/(ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.mean())
print(" \n\n\n the normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData)
print(" \nthe variance of normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.var())
print(" the stand diviation of normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData.std())

##################################################################################

print(" \n\n\n the mean of Thin-Wall Energy Barrier(experiemt) is : ", ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.mean())
print("  the variance of Thin-Wall Energy Barrier(experiemt) is : ", ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.var())
print("  the stand diviation of Thin-Wall Energy Barrier(experiemt) is : ", ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.std())

normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData = ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData/(ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.mean())
print(" \n\n\n the normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)
print(" \nthe variance of normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.var())
print(" the stand diviation of normalized Thin-Wall Energy Barrier(experiemt) is : ", normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData.std())

    
##############################################################################
# Scatter plot Energy Barriers and Ratio of the Lotynk's Experiment parameters
##############################################################################

n = np.arange(0, lengthPressure, 1)
plot1.scatter(pressure, normalized_ratio_fAfAB_for_LotynkExperimentData)
plot1.scatter(pressure, normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)
plot1.xlabel(r'$p/bar$');
plot1.legend([r"$\frac{|f_{A}|}{f_{AB}}, StandDiviation = 0.1185 $" , r"$\frac{E^{Thin-Wall}_{c}}{{K_b}T}, \sigma_{AB} = 0.71|f_{B}|{\xi_{GL}}, , StandDiviation = 0.3037$"])
plot1.title(" Normalized f_A/f_AB and Thin-Wall E_c, Data Distributions and Stand Diviations ")
plot1.savefig('normalized_fAfAB_ThinWall_EnergyBarrier_with_Lotynk_Experiment.pdf')
plot1.show()


plot1.clf()
plot1.cla()
plot1.close()

n = np.arange(0, lengthPressure, 1)
plot1.scatter(pressure, ratio_fAfAB_for_LotynkExperimentData)
plot1.xlabel(r'$p/bar$');plot1.ylabel(r'$|f_{A}|/|f_{A}-f_{B}|$')
plot1.legend([r"$\frac{|f_{A}|}{|f_{A}-f_{B}|}, mean = 302.9, std = 36.0 $"])
plot1.title(" AB nucleation data from Lotnyk et al 2021 ")
plot1.savefig('fAfAB_Lotynk_Experiment.pdf')
plot1.show()

plot1.clf()
plot1.cla()
plot1.close()             


n = np.arange(0, lengthPressure, 1)
plot1.scatter(pressure, ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)
plot1.xlabel(r'$p/bar$');
plot1.legend([r"$\frac{E^{Thin-Wall}_{c}}{{K_b}T}, MeanValue =88078802, StandDiviation = 26767375$"])
plot1.title(" Thin-wall energy barrier Data Distributions and Stand Diviations ")
plot1.savefig('Thin_wall_energy_barrier_Lotynk_Experiment.pdf')
plot1.show()

plot1.clf()
plot1.cla()
plot1.close()             

             
plot1.scatter(n, normalized_ratio_fAfAB_for_LotynkExperimentData)
plot1.scatter(n, normalized_ratio_ExperimentTension_for_LotynkExperimentData)
plot1.scatter(n, normalized_ThinWall_Evaluation_of_EnergyBarrier_ExperimentTension_LotynkExperimentData)
plot1.scatter(n, normalized_ThinWall_Evaluation_of_EnergyBarrier_fAfAB_LotynkExperimentData)
plot1.legend([r"$\frac{|f_{A}|}{f_{AB}}$" , r"$\frac{0.71 |f_{B}|}{f_{AB}}$", r"$E^{Thin-Wall}_{c}, \sigma_{AB} = 0.71|f_{B}|{\xi_{GL}}$", r"$E^{Thin-Wall}_{c}, \sigma_{AB} = |f_{A}|{\xi_{GL}}$"])
plot1.show()


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



# ###################################################
# # contour plot Rc of the |f_A| xi evaluation
# ###################################################

# # LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3, 7*10**3, 9*10**3, 10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8]
# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3] 
# X, Y = np.meshgrid(Temperature, pressure)
# # fig, ax = plot1.subplots()
# cs1 = plot1.contourf(X*(10**3), Y, ThinWall_Estimate_Rc_fAxi, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
# cb = plot1.colorbar(cs1)
# cb.set_ticks([10**0, 5*10**0, 10**1, 5*10**1, 10**2, 5*10**2, 10**3])

# #cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi} /\xi_{GL}$')
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

# # plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_GL}})^{3}}{{{\Delta}f_{AB}}^{2}} $')
# plot1.title(r'$R^{thin-wall}_c/{\xi_{GL}} \sim \frac{2 |f_A|}{{\Delta}f_{AB}} /\xi_{GL}$')
# plot1.savefig('Distribution_of_Rc_ThinWall_Evaluation_Contour_fAxi_SCModule_versionI.pdf');

# plot1.show()

# ####################################################################
# # contour plot of Rc by 0.7 |f_B|xi tension
# ####################################################################


# X, Y = np.meshgrid(Temperature, pressure)
# # fig, ax = plot1.subplots()
# cs1 = plot1.contourf(X*(10**3), Y, ThinWall_Estimate_Rc_ExpermentTension, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
# cb = plot1.colorbar(cs1)
# cb.set_ticks([10**0, 5*10**0, 10**1, 5*10**1, 10**2, 5*10**2, 10**3])
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

# plot1.title(r'$R^{thin-wall}_{c} \sim \frac{2 \times 0.71 |f_B|}{{\Delta}f_{AB}} /\xi_{GL}$')
# plot1.savefig('Distribution_of_Rc_ThinWall_Evaluation_ExperimentalTension_Contour_SCModule_versionI.pdf');

# plot1.show()



# #################################################################
# # contour plot of temperature-pressure dependent coherent length
# #################################################################




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



