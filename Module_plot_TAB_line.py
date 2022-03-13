
# This script is used for calculating the Thin-wall evaluation of radius and energy barrier of critical bubble.
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


import Module_SCCO_V02 as SC # strong coupling correction module
import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math
from matplotlib import ticker, cm
import csv

# main code starts from here

m = 1;s = 1; J = 1; Kelvin = 1; bar = 1# Length unit, Time unit, Energy unit, Temperature unit, pressure unit

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
# calculate free energy under different temperature     

stepT = 0.01*(10**-3) 
Temperature = np.arange(1.5*(10**-3), 2.40*(10**-3), stepT) #Kelvin

stepPressure = 0.01*bar
# pressure = np.arange(20.0, 33.9*bar, stepPressure)
pressure = np.arange(19.9, 33.99*bar+stepPressure, stepPressure)

X, Y = np.meshgrid(Temperature, pressure)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)

# lots of saving matrices

TempretureDependent_GL_CoherentLength = np.zeros((lengthPressure,lengthT)) # save the temperature dependent coherent length

EnergyDensity_Difference_fABGL = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fAB

EnergyDensity_fA = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fA

EnergyDensity_fB = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fB


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
    
    for iT in range(0, lengthT, 1):
        #indexDelta = math.floor(delta/stepDelta)
        indexT = iT
        print('indexT is',indexT, ' Temperature is ', Temperature[indexT])
       
        t = Temperature[indexT]/Tcp
        print('temperatureis:, ',t)

        xiGLWeakCoupling = math.sqrt((7.0*zeta3)/20.0)*xi0p # weak coupling GL coherent length in PRB 101 024517

         

        if t >= 1:

           print(" bro, we just got temperature at Tc, save a np.nan. ")
      
          
           # save the energy density difference
           EnergyDensity_Difference_fABGL[indexP,indexT] = np.nan

           # save the energy density difference
           EnergyDensity_fA[indexP,indexT] = np.nan
           # save the energy density difference
           EnergyDensity_fB[indexP,indexT] = np.nan

           # masked GL coherent length
           TempretureDependent_GL_CoherentLength[indexP,indexT] = np.nan
           

           
        else:

           xitp = xiGLWeakCoupling/math.sqrt(1-t) # temperature dependent coherent length  

           # save the temperature-dependent coherent length
           TempretureDependent_GL_CoherentLength[indexP,indexT] = xitp
           
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

           # save the energy density difference
           EnergyDensity_Difference_fABGL[indexP,indexT] = DiffFABGL

           # save the energy density difference
           EnergyDensity_fA[indexP,indexT] = fAGL
           # save the energy density difference
           EnergyDensity_fB[indexP,indexT] = fBGL
           
           
           print('fAGL in SI unit is:', fAGL)
           print('fBGL in SI unit is:', fBGL)
           print('DiffFABGL in SI unit is:', DiffFABGL)   

    

# Plot the TAB line
# cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

# plot1.show()




