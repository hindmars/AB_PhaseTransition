
# This script comes from the Version I, which is for calculating the Thin-wall evaluation of radius and energy barrier of critical bubble.
# The basic idea is using A. J. Leggtt's simply considerations when f_A - f_B is small and
# R_c / \xi, R_c/t are both very tiny, and close to T_AB
# For details, checking his Journal of Low Temperature vol.87 571 1992 and Yip & Leggtt review
# in Helium Three 1990.

# the critical radius are plotted as countours, and compared with the experiment data from Lotynk's experiment.

# Mark suggested another evaluation about R_c in 19th Nov 2021, see the email in details
# the basic idea is R_c ~ ((f_A \xi_GL)/ \Deltaf_AB)/\xi_GL = f_A/ \Deltaf_AB
# temperature region from 0.0 mk to 2.4 mK

# Now is 29th Nov. 2021, the accidentaly coincide of f_A/ f_AB seems is kind of clue of " thick wall " bubble
# from Enqvist paper (PRD 1992 with Kari Lummukainen). Mark and me had a talk in this afternoon for this.
# If there is E_x/KT ~ |f_A|/f_AB and \Gamma ~ exp^(-E_x/KT), then Lotynk's expermental data will drop in to same
# contour of |f_A|/f_AB as long as E_x/KT(|f_A/f_AB|) is monitonic.

# So, this Version II script focuses on plot contours for |f_A|/f_AB, 2 |f_A|/f_AB, etc.
# And because the numpy.array seems induced some calculation error, which seems it can't directly
# take the elements lever calculation, I will still us loop for every elements ot make sure results are relaible.

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
# calculate free energy under different temperature     

stepT = 0.001*(10**-3) 
Temperature = np.arange(2.0*(10**-3), 2.40*(10**-3), stepT) #Kelvin

stepPressure = 0.001*bar
# pressure = np.arange(20.0, 33.9*bar, stepPressure)
pressure = np.arange(19.9, 22.99*bar+stepPressure, stepPressure)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)

# lots of saving matrices
DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the energy differene, which be scaled by |fBGL|

ThinWall_Evaluation_of_EnergyBarrier_DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of energy barrier with experiment tension \sigma = 0.7 \xi f_{B}

ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of energy barrier with insteading tension

ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of energy barrier with insteading tension by |f_A| \xi




TempretureDependent_GL_CoherentLength = np.zeros((lengthPressure,lengthT)) # save the temperature dependent coherent length

EnergyDensity_Difference_fABGL = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fAB

EnergyDensity_fA = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fA

EnergyDensity_fB = np.zeros((lengthPressure,lengthT)) # save the SI unit \Delta fB

ThinWall_Estimate_Rc_fAxi = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of energy barrier with insteading tension by |f_A| \xi

ThinWall_Estimate_Rc_ExpermentTension = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of R_{c} with experiment tension \sigma = 0.7 \xi f_{B}



ratio_fA_fAB_PreIndexIs1 = np.zeros((lengthPressure,lengthT)) # save the |f_A|/f_AB

ratio_fA_fAB_PreIndexIs2 = np.zeros((lengthPressure,lengthT)) # save the 2 |f_A|/f_AB

ratio_fA_fAB_PreIndexIs3 = np.zeros((lengthPressure,lengthT)) # save the 3 |f_A|/f_AB

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
      
           ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi[indexP,indexT] = np.nan
           ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension[indexP,indexT] = np.nan

           ThinWall_Estimate_Rc_fAxi[indexP,indexT] = np.nan
           ThinWall_Estimate_Rc_ExpermentTension[indexP,indexT] = np.nan


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

           
           ################################################################
           # calculation for thin-wall evaluation of energy barrier
           ################################################################
           numerator1 = (abs(fAGL)*xitp)**3 # tension for |fA| xi evaluation
           numerator2 = (0.71*abs(fBGL)*xitp)**3 # tension from experiment closed TAB
           dinorminator = DiffFABGL**2 
           KbT = kb*Temperature[indexT]
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
           # calculation for |f_A|/f_AB, 2 |f_A|/f_AB, 3 |f_A|/f_AB
           ###################################################################
           
           fAfABPre1 = (1.0*abs(fAGL)*1)/DiffFABGL

           fAfABPre2 = (2.0*abs(fAGL)*1)/DiffFABGL

           fAfABPre3 = (3.0*abs(fAGL)*1)/DiffFABGL

        
           
           if DiffFABGL > 0:

              
              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi[indexP,indexT] = EnergyBarrier_fAxi_KbT # in unit of |f_B| \xi_GL^3

              print(" thin wall energy barrier in KbT is ", ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi[indexP,indexT])

              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension[indexP,indexT] = EnergyBarrier_experimentTension_KbT # in unit of |f_B| \xi_GL^3

              print(" thin wall energy barrier in KbT is ", ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension[indexP,indexT])

              
              

              ratio_fA_fAB_PreIndexIs1[indexP,indexT] = fAfABPre1

              print(" 1*|f_A|/f_AB is ", ratio_fA_fAB_PreIndexIs1[indexP,indexT])

              ratio_fA_fAB_PreIndexIs2[indexP,indexT] = fAfABPre2

              print("  2*|f_A|/f_AB is ", ratio_fA_fAB_PreIndexIs2[indexP,indexT])

              ratio_fA_fAB_PreIndexIs3[indexP,indexT] = fAfABPre3

              print("  3*|f_A|/f_AB is ", ratio_fA_fAB_PreIndexIs3[indexP,indexT])


           else:

              print(" bro, you get negative logarithm, save a np.nan. ")

              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi[indexP,indexT] = np.nan

              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension[indexP,indexT] = np.nan

              ThinWall_Estimate_Rc_fAxi[indexP,indexT] = np.nan

              ThinWall_Estimate_Rc_ExpermentTension[indexP,indexT] = np.nan



              ratio_fA_fAB_PreIndexIs1[indexP,indexT] = np.nan

              ratio_fA_fAB_PreIndexIs2[indexP,indexT] = np.nan

              ratio_fA_fAB_PreIndexIs3[indexP,indexT] = np.nan

            
        

    print('Thin-wall evaluation of Energy barrier by |fA|xiin KbT is: ', ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_fAxi[indexP,:])
    print('Thin-Wall evaluation of Energy barrier by 0.7 |fB|xi in KbT is: ', ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL_ExperimentTension[indexP,:])

    print('Thin-wall evaluation of Rc by |fA|xiin KbT is: ', ThinWall_Estimate_Rc_fAxi[indexP,:])
    print('Thin-Wall evaluation of Rc by 0.7 |fB|xi in KbT is: ', ThinWall_Estimate_Rc_ExpermentTension[indexP,:])

    print('the GL coherent length in this pressure is: ', TempretureDependent_GL_CoherentLength[indexP,:]*(10**9))
    


###################################################
# contour plot |f_A|/f_AB 
###################################################

# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3, 7*10**3, 9*10**3, 10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8]
# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3]
LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3] 
X, Y = np.meshgrid(Temperature, pressure)
# fig, ax = plot1.subplots()
cs1 = plot1.contourf(X*(10**3), Y, ratio_fA_fAB_PreIndexIs1, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
cb = plot1.colorbar(cs1)
cb.set_ticks([10**0, 2*10**0, 3*10**0, 4*10**0, 5*10**0, 6*10**0, 7*10**0, 8*10**0, 9*10**0, 10**1, 2*10**1, 3*10**1, 4*10**1, 5*10**1, 6*10**1, 7*10**1, 8*10**1, 9*10**1, 10**2, 2*10**2, 3*10**2, 4*10**2, 5*10**2, 6*10**2, 7*10**2, 8*10**2, 9*10**2, 10**3, 2*10**3, 3*10**3, 4*10**3, 5*10**3])

#cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi} /\xi_{GL}$')
# plot1.legend(loc='lower right')

# Plot the TAB line
cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

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

plot1.scatter(fuckData1, fuckData2)
#matplotlib.pyplot.scatter(fuckData1, fuckData2)

# plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_GL}})^{3}}{{{\Delta}f_{AB}}^{2}} $')
plot1.title(r'$\frac{|f_A|}{{\Delta}f_{AB}}$')
plot1.savefig('Distribution_of_fA_fAB_PreIndex1_SCModule.pdf');

plot1.show()

###################################################
# contour plot 2 |f_A|/f_AB 
###################################################

# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3, 7*10**3, 9*10**3, 10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8]
# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3]
LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3] 
X, Y = np.meshgrid(Temperature, pressure)
# fig, ax = plot1.subplots()
cs1 = plot1.contourf(X*(10**3), Y, ratio_fA_fAB_PreIndexIs2, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
cb = plot1.colorbar(cs1)
cb.set_ticks([10**0, 2*10**0, 3*10**0, 4*10**0, 5*10**0, 6*10**0, 7*10**0, 8*10**0, 9*10**0, 10**1, 2*10**1, 3*10**1, 4*10**1, 5*10**1, 6*10**1, 7*10**1, 8*10**1, 9*10**1, 10**2, 2*10**2, 3*10**2, 4*10**2, 5*10**2, 6*10**2, 7*10**2, 8*10**2, 9*10**2, 10**3, 2*10**3, 3*10**3, 4*10**3, 5*10**3])

#cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi} /\xi_{GL}$')
# plot1.legend(loc='lower right')

# Plot the TAB line
cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

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

plot1.scatter(fuckData1, fuckData2)
#matplotlib.pyplot.scatter(fuckData1, fuckData2)

# plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_GL}})^{3}}{{{\Delta}f_{AB}}^{2}} $')
plot1.title(r'$\frac{2 |f_A|}{{\Delta}f_{AB}}$')
plot1.savefig('Distribution_of_fA_fAB_PreIndex2_SCModule.pdf');

plot1.show()


###################################################
# contour plot 3 |f_A|/f_AB 
###################################################

# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3, 7*10**3, 9*10**3, 10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6, 7*10**6, 9*10**6, 10**7, 3*10**7, 5*10**7, 7*10**7, 9*10**7, 10**8]
# LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3]
LLLLL = [10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3] 
X, Y = np.meshgrid(Temperature, pressure)
# fig, ax = plot1.subplots()
cs1 = plot1.contourf(X*(10**3), Y, ratio_fA_fAB_PreIndexIs3, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
cb = plot1.colorbar(cs1)
cb.set_ticks( [10**0, 2*10**0, 3*10**0, 4*10**0, 5*10**0, 6*10**0, 7*10**0, 8*10**0, 9*10**0, 10**1, 2*10**1, 3*10**1, 4*10**1, 5*10**1, 6*10**1, 7*10**1, 8*10**1, 9*10**1, 10**2, 2*10**2, 3*10**2, 4*10**2, 5*10**2, 6*10**2, 7*10**2, 8*10**2, 9*10**2, 10**3, 2*10**3, 3*10**3, 4*10**3, 5*10**3])

#cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi} /\xi_{GL}$')
# plot1.legend(loc='lower right')

# Plot the TAB line
cs2 = plot1.contour(X*(10**3), Y, EnergyDensity_Difference_fABGL, levels=[0.0], colors='red'); # plot1.clabel(Cs2, inline=True, fontsize=8.5, colors='r')

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

plot1.scatter(fuckData1, fuckData2)
#matplotlib.pyplot.scatter(fuckData1, fuckData2)

# plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_GL}})^{3}}{{{\Delta}f_{AB}}^{2}} $')
plot1.title(r'$\frac{3 |f_A|}{{\Delta}f_{AB}}$')
plot1.savefig('Distribution_of_fA_fAB_PreIndex3_SCModule.pdf');

plot1.show()
