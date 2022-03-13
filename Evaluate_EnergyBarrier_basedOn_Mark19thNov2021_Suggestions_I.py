
# This script is used for calculating the Thin-wall evaluation of radius of critical bubble
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

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version.2 of ebergy difference code, in which the pico-Joule(pJ 10^-12) unit is used 

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569;



import matplotlib.pyplot as plot1
import matplotlib
import numpy as np

from math import pi
import math
from numpy import ma
from matplotlib import ticker, cm
import csv


def region_judgement(p):
    if (p >= 0) and (p < 2):
         pk=0;pk1=2;k=0
    if (p >= 2) and (p < 4):
         pk=2;pk1=4;k=1
    if (p >= 4) and (p < 6):
         pk=4;pk1=6;k=2
    if (p >= 6) and (p < 8):
         pk=6;pk1=8;k=3
    if (p >= 8) and (p < 10):
         pk=8;pk1=10;k=4
    if (p >= 10) and (p < 12):
         pk=10;pk1=12;k=5
    if (p >= 12) and (p < 14):
         pk=12;pk1=14;k=6
    if (p >= 14) and (p < 16):
         pk=14;pk1=16;k=7
    if (p >= 16) and (p < 18):
         pk=16;pk1=18;k=8     
    if (p >= 18) and (p < 20):
         pk=18;pk1=20;k=9
    if (p >= 20) and (p < 22):
         pk=20;pk1=22;k=10
    if (p >= 22) and (p < 24):
         pk=22;pk1=24;k=11
    if (p >= 24) and (p < 26):
         pk=24;pk1=26;k=12
    if (p >= 26) and (p < 28):
         pk=26;pk1=28;k=13
    if (p >= 28) and (p < 30):
         pk=28;pk1=30;k=14     
    if (p >= 30) and (p < 32):
         pk=30;pk1=32;k=15
    if (p >= 32) and (p <= 34):
         pk=32;pk1=34;k=16

    return [pk, pk1, k] #return list


def intepolation_presure_SC(pk,pk1,p,fk,fk1):
    K = (fk1-fk)/(pk1-pk);print("K =",K, "pk1 =",pk1, "pk =",pk, "fk1 =",fk1, "fk =",fk)
    fip = K * (p-pk) + fk;print("fip =",fip)

    # return the cip or Tcp
    return fip
         
class BETA:

     # class of beta object
    

     def __init__(self,name):
         self.name = name
         print(" beta oject is crated ")

     def c1_function(self,pressure,pk,pk1,k):
         c1 = [-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413]
         self.c1p =  intepolation_presure_SC(pk,pk1,pressure,c1[k],c1[k+1])

     def c2_function(self,pressure,pk,pk1,k):
         c2 = [-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645]
         self.c2p =  intepolation_presure_SC(pk,pk1,pressure,c2[k],c2[k+1])

     def c3_function(self,pressure,pk,pk1,k):
         c3 = [-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268]
         self.c3p =  intepolation_presure_SC(pk,pk1,pressure,c3[k],c3[k+1])

     def c4_function(self,pressure,pk,pk1,k):
         c4 = [-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518]
         self.c4p =  intepolation_presure_SC(pk,pk1,pressure,c4[k],c4[k+1])

     def c5_function(self,pressure,pk,pk1,k):
         c5 = [-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815]
         self.c5p =  intepolation_presure_SC(pk,pk1,pressure,c5[k],c5[k+1])

     def tc_function(self,pressure,pk,pk1,k):
         Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]
         self.tcp =  intepolation_presure_SC(pk,pk1,pressure,Tc[k],Tc[k+1])

     def mStar_function(self,pressure,pk,pk1,k):
         Ms = [2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82] # in unit of helium-3 atom
         self.ms =  intepolation_presure_SC(pk,pk1,pressure,Ms[k],Ms[k+1])

     def vFermi_function(self,pressure,pk,pk1,k):
         VF = [59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23]
         self.vf = intepolation_presure_SC(pk,pk1,pressure,VF[k],VF[k+1])
     





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

BetaObject = BETA('betaAndTc')

#p = 21.1 # pressure, bar
# pressureArray = np.arange(0.0, 34.0, 2.0)

#for p in pressureArray:
    
###################################################
# calculate free energy under different temperature     

stepT = 0.005*(10**-3) 
Temperature = np.arange(2.0*(10**-3), 2.40*(10**-3), stepT) #Kelvin

stepPressure = 0.005*bar
# pressure = np.arange(20.0, 33.9*bar, stepPressure)
pressure = np.arange(19.9, 22.99*bar+stepPressure, stepPressure)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)

fBGL_array = np.zeros((lengthPressure,lengthT))
fAGL_array = np.zeros((lengthPressure,lengthT))
DiffFABGL = np.zeros((lengthPressure,lengthT)) # save data of the energy difference

DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the energy differene, which be scaled by |fBGL|

ThinWall_Estimate_Rc_DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of R_{c} with experiment tension \sigma = 0.7 \xi f_{B}

ThinWall_Evaluation_of_EnergyBarrier_DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the thin-wall evaluation of R_{c} with experiment tension \sigma = 0.7 \xi f_{B}

ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL = np.zeros((lengthPressure,lengthT))

LeggttEstimate_Rc_DiffFABGLScaled_logarithim = np.zeros((lengthPressure,lengthT)) # save logarithm data of the thin-wall evaluation of R_{c} with experiment tension \sigma = 0.7 \xi f_{B}

for iP in range(0, lengthPressure, 1):
    print('\n\n Now P is:', pressure[iP], '\n\n')
    #indexT = math.floor(T/stepT)
    indexP = iP
    print('indexP is ',indexP)

    p = pressure[iP]
    judgementList = region_judgement(p)

    print('judgementList is', judgementList)

    print('low pressure is: ',judgementList[0], ',high pressure is: ',judgementList[1], ',interpolation region is: ', judgementList[2])

    pk = judgementList[0]; pk1 = judgementList[1]; k = judgementList[2]

    BetaObject.c1_function(p,pk,pk1,k);c1p = BetaObject.c1p
    BetaObject.c2_function(p,pk,pk1,k);c2p = BetaObject.c2p
    BetaObject.c3_function(p,pk,pk1,k);c3p = BetaObject.c3p
    BetaObject.c4_function(p,pk,pk1,k);c4p = BetaObject.c4p
    BetaObject.c5_function(p,pk,pk1,k);c5p = BetaObject.c5p
    BetaObject.tc_function(p,pk,pk1,k);Tcp = (BetaObject.tcp)*(10**(-3))
    BetaObject.mStar_function(p,pk,pk1,k); mEffective = (BetaObject.ms)*m3
    BetaObject.vFermi_function(p,pk,pk1,k); vFermi = BetaObject.vf
    

    c245p = c2p + c4p + c5p;c12p = c1p + c2p;c345p = c3p + c4p + c5p

    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,'\n\n')

    N0 = ((mEffective**(2))*vFermi)/((2*pi*pi)*(hbar**(3))) # energy density of Fermi surface

    print('\npressure is, ',p,' effective mass is, ', mEffective, ' Fermi velocity is,', vFermi, ' N(0) is ',N0,'\n\n')
    
    for iT in range(0, lengthT, 1):
        #indexDelta = math.floor(delta/stepDelta)
        indexT = iT
        print('indexT is',indexT)
       
        t = Temperature[indexT]/Tcp
        print('temperatureis:, ',t)

        if t > 1:

           print(" bro, we just got temperature at Tc, save a np.nan. ")

           DiffFABGL[indexP,indexT] = np.nan
           DiffFABGLScaled[indexP,indexT] = np.nan

           fBGL_array[indexP,indexT] = fBGLRed
           fAGL_array[indexP,indexT] = fAGLRed

           
        else:    

           alphaRed = (1/3)*(t-1)

           beta245Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(2+t*c245p) # A Phase
           beta12345Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(3.0*(1+t*c12p) + (2+t*c345p)) # B phase

           deltaA = math.sqrt((-alphaRed)/(4*beta245Red)) # A -> B

           deltaB = math.sqrt((-alphaRed)/(2*beta12345Red)) # A -> B
        
           fAGLRed = alphaRed*2*(deltaA**2) + 4*beta245Red*(deltaA**4)
 
#       fBGLRed = alphaRed*3*(delta**2) + 3*beta12345Red*(delta**4)
           fBGLRed = alphaRed*3*(deltaB**2) + 3*beta12345Red*(deltaB**4)

           fBGL = fBGLRed * N0
           print('fBGL is:', fBGL)
           absfBGLxiGL3 = abs(fBGL) * ((20*(10**(-9)))**3) # 20 nm (10^⁻9) coherent around 20 bar
           print('|fBGL| is:', absfBGLxiGL3)
           absfBGLxiGL3_KbT = absfBGLxiGL3/(kb*t*Tcp)
           print('absfBGLxiGL3_KbT is:', absfBGLxiGL3_KbT)


           
           DiffFABGLScaled[indexP,indexT] = (fAGLRed - fBGLRed)/abs(fBGLRed)
           TempData1 = (fAGLRed - fBGLRed)/abs(fBGLRed) # energy difference in the unit of f_B
           # TempData2 = abs(fAGLRed)/abs(fBGLRed) # geberate f_A with unit of |f_B|
           TempData2 = 0.7 # numerical value of surface
           
           if TempData1 > 0:

              
              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGLScaled[indexP,indexT] =(16/3)*pi*((TempData2**3)/TempData1) # in unit of |f_B| \xi_GL^3
              ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL[indexP,indexT] =((16/3)*pi*((TempData2**3)/(TempData1)**2)) * absfBGLxiGL3_KbT # in unit of |f_B| \xi_GL^3
              # LeggttEstimate_Rc_DiffFABGLScaled_logarithim[indexP,indexT] = math.log((2*0.7)/TempData)
              print(" thin wall Rc is ", ThinWall_Estimate_Rc_DiffFABGLScaled[indexP,indexT])

           else:

              print(" bro, you get negative logarithm, save a np.nan. ")

              ThinWall_Estimate_Rc_DiffFABGLScaled[indexP,indexT] = np.nan
              #vLeggttEstimate_Rc_DiffFABGLScaled_logarithim[indexP,indexT] = np.nan  
                
          # DiffFABGL[indexP,indexT] = DiffFABGLScaled[indexP,indexT]*abs(fBGLRed*N0)
        
           fBGL_array[indexP,indexT] = fBGLRed
           fAGL_array[indexP,indexT] = fAGLRed

    print('fBGLRed_Delta is:', fBGL_array[indexP,:])
    print('fAGLRed_Delta is:', fAGL_array[indexP,:])
    print('difference of fGLRed_Delta is:', DiffFABGL[indexP,:])
    print('difference between A & B fGLRed_Delta is:', DiffFABGLScaled[indexP,:])
    print('Mark suggested evaluation of Rc is: ', ThinWall_Estimate_Rc_DiffFABGLScaled[indexP,:])
#    plo
#    plot1.plot(Delta, fGL_array[indexT,:], 'o-'); plot1.ylabel('ev^2'); plot1.xlabel('ev')
#    plot1.show()
    

# # Plot Free energy difference
# for iP in range(0, lengthPressure, 1):
#     #indexT = math.floor(T/stepT)
#     indexP = iP
#     plot1.plot(Temperature*1000, DiffFABGL[indexP,:], '-.'); plot1.ylabel(r'$(f{A}_{GL}-f^{B}_{GL}).N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')
    
# plot1.show()    
# plot1.clf()
# plot1.cla()
# plot1.close()    

# # Plot Free Energies for A and B 
# for iP in range(0, lengthPressure, 1):
#     #indexT = math.floor(T/stepT)
#     indexP = iP
   
#     plot1.plot(Temperature*1000, fAGL_array[indexP,:], '--'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')
#     plot1.plot(Temperature*1000, fBGL_array[indexP,:], '-'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')

# plot1.show()    
 
# plot1.clf()
# plot1.cla()
# plot1.close()    


# density and contour plot of (fAGL - fBGL)/|fBGL|
# Levels = np.arange(0.0, 0.2, 0.01)
# DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, DiffFABGLScaled);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)
# Cplot = plot1.contour(Temperature*1000, pressure, DiffFABGLScaled, levels=Levels, colors='black');plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
# Cplot.collections[0].set_label(r'$(f_{A}-f_{B})/|f_{B}|$')
# plot1.legend(loc='lower right')
# plot1.savefig('DensityPlot_FreeEnergyDiff_Scaled_SI_unit.pdf');
# plot1.show()

# plot1.clf()
# plot1.cla()
# plot1.close()

# calculate Leggtt Estimate about Rc
Levels = [0.0]
# DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, LeggttEstimate_Rc_DiffFABGLScaled);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)
# Cplot = plot1.contour(Temperature*1000, pressure, LeggttEstimate_Rc_DiffFABGLScaled, levels=Levels, colors='black');plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
Cplot = plot1.contour(Temperature*1000, pressure, ThinWall_Evaluation_of_EnergyBarrier_DiffFABGLScaled, colors='black');plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
# Cplot.collections[0].set_label(r'$2 {\sigma^{experiment}_{AB}}/{\Delta}f_{AB}$')
# plot1.legend(loc='lower right')
# plot1.savefig('Distribution_of_Rc_Leggtt_Estimateion.pdf');
plot1.show()

plot1.clf()
plot1.cla()
plot1.close()

# calculate Leggtt Estimate about Rc, Logarithm corrdinate
Levels = np.arange(1, 80, 1)
bounds = [1*10**0, 2*10**0, 3*10**0, 4*10**0, 5*10**0] 
DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, ThinWall_Evaluation_of_EnergyBarrier_DiffFABGLScaled);
# plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot, extend='max', ticks=bounds)
#DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, LeggttEstimate_Rc_DiffFABGLScaled,
#                               norm=matplotlib.colors.LogNorm(vmin=LeggttEstimate_Rc_DiffFABGLScaled.min(), vmax=LeggttEstimate_Rc_DiffFABGLScaled.max()));
# plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot, extend='max')
plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)


# # Cplot = plot1.contour(Temperature*1000, pressure, LeggttEstimate_Rc_DiffFABGLScaled_logarithim, levels=Levels, colors='black');plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
# Cplot.collections[0].set_label(r'$2 {\sigma^{experiment}_{AB}}/{\Delta}f_{AB}$')
# plot1.legend(loc='lower right')
# plot1.savefig('Distribution_of_Rc_Leggtt_Estimateion_LogarithmCoordinate_Logarithim.pdf');
plot1.show()

plot1.clf()
plot1.cla()
plot1.close()

# calculate Leggtt Estimate about Rc, Logarithm corrdinate
# Levels = np.arange(, , 1)
# Cplot = plot1.contourf(Temperature*1000, pressure, LeggttEstimate_Rc_DiffFABGLScaled, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
# Cplot.collections[0].set_label(r'$2 {\sigma^{experiment}_{AB}}/{\Delta}f_{AB}$')
# plot1.legend(loc='lower right')
# plot1.savefig('Distribution_of_Rc_Leggtt_Estimateion_LogarithmCoordinate_Logarithim_Contour.pdf');
# plot1.show()

# plot1.clf()
# plot1.cla()
# plot1.close()

#LLLLL = [1*10**0, 3*10**0, 5*10**0, 7*10**0, 9*10**0, 10**1, 3*10**1, 5*10**1, 7*10**1, 9*10**1, 10**2, 3*10**2, 5*10**2, 7*10**2, 9*10**2, 10**3, 3*10**3, 5*10**3, 7*10**3, 9*10**3, 10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5]

LLLLL = [10**4, 3*10**4, 5*10**4, 7*10**4, 9*10**4, 10**5, 3*10**5, 5*10**5, 7*10**5, 9*10**5, 10**6, 3*10**6, 5*10**6] 
X, Y = np.meshgrid(Temperature, pressure)
# fig, ax = plot1.subplots()
cs = plot1.contourf(X*(10**3), Y, ThinWall_Evaluation_of_EnergyBarrier_DiffFABGL, locator=ticker.LogLocator(), cmap=cm.PuBu_r, levels=LLLLL);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');
#cs.collections[0].set_label(r'${R^{thin-wall}_{c}}/{\xi}=2{\sigma^{experiment}_{AB}}/{\Delta}f_{AB}{\xi}$')
# plot1.legend(loc='lower right')

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

# #atplotlib.pyplot.scatter(list(list1[0]),list(list1[1]))
# #matplotlib.pyplot.show()

plot1.scatter(fuckData1, fuckData2)
#matplotlib.pyplot.scatter(fuckData1, fuckData2)

plot1.title(r'${\Delta}G^{thin-wall}_{Barrier}/(k_{B}T) \sim \frac{16{\pi}}{3} \frac{({0.7|f_B|{\xi_GL}})^{3}}{{{\Delta}f_{AB}}^{2}} $')
plot1.colorbar(cs)
plot1.savefig('Distribution_of_EnergyBarrier_ThinWall_Evaluation_Logarithim_Contour_versionIII.pdf');
# plot1.show()

plot1.show()



# for iP in [70, 80, 90, 100, 110, 120]:
#     print(" pressure is: ",pressure[iP])
#     #indexT = math.floor(T/stepT)
#     indexP = iP
   
#     plot1.plot(Temperature*1000, DiffFABGLScaled[indexP,:], '-.'); plot1.ylabel(r'$(f{A}_{GL}-f^{B}_{GL})$/$|f^{B}|$'); plot1.xlabel(r'$T$/mK')

# plot1.savefig('FreeEnergyDiff_Scaled.pdf');
# plot1.show()

