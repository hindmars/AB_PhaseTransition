# this script is for plot-check the GL free energy density under given pressure
# the strong coupling corrected \beta coefficients are taken from PRB. 92. 144515, 

# the out put plot coincide with result in Fig.1 of PRB. 92. 144515 



# zeta3 = 1.2020569;



import matplotlib.pyplot as plot1
import numpy as np

from math import pi
import math


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
    if (p >= 32) and (p < 34):
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





# main code starts from here

# pressure = np.arange(0.0, 34.05, 0.05);

m = 1;s = 1; eV = 1; Kelvin = 1 # Length unit, Time unit, Energy unit
zeta3 = 1.2020569;
kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1
c = 2.99792458*(10**(8)) # speed of light, m.s^-1
hbar = 6.582119569*(10**(-16)) #plank constant, eV.s
u = 9.3149410242*(10**(8))*eV*(c**(-2)) #atomic mass unit, Dalton, eV.c^-2
m3 = 3.016293*u #mass of helium3 atom
ms = 5.02*m3 #effective mass under 24bar
vf = 36.53 #fermi velosity under 24bar, m.s^1

kf = (vf*ms)/hbar
Nf = (ms*kf)/(pi*pi*hbar*hbar)

C1=[];C2=[];C3=[];C4=[];C5=[];Tc=[];


###########################################
# build object, and check the intepolatinon 

BetaObject = BETA('betaAndTc')

p = 21.1 # pressure, bar
# pressureArray = np.arange(0.0, 34.0, 2.0)

#for p in pressureArray:
    
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

c245p = c2p + c4p + c5p

c12p = c1p + c2p

c345p = c3p + c4p + c5p

print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,'\n\n\n')



###################################################
# calculate free energy under different temperature     

stepT = 0.5*(10**-3) 
Temperature = np.arange(0.0, 2.5*(10**-3), stepT) #Kelvin
#Temperature = [1.8*(10**(-3)), 2.275*(10**(-3))]

stepDelta = 0.1*kb*Tcp
Delta = np.arange(0.0, 5*kb*Tcp, stepDelta)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Delta is',Delta,'\n length of Delta is ', len(Delta))
lengthDelta = len(Delta)

fBGL_array = np.zeros((lengthT,lengthDelta))
fAGL_array = np.zeros((lengthT,lengthDelta))

for iT in range(0, lengthT, 1):
    print('\n\n Now T is:', Temperature[iT], '\n\n')
    #indexT = math.floor(T/stepT)
    indexT = iT
    print('indexT is ',indexT)
    
    for iD in range(0, lengthDelta, 1):
        #indexDelta = math.floor(delta/stepDelta)
        indexDelta = iD
        print('indexDelta is',indexDelta)
       
        t = Temperature[indexT]/Tcp
        #print('pressure is:, ',p, 'Tcp is:,', tcp)

        delta = Delta[indexDelta]
          
        alphaRed = (1/3)*(t-1)
        beta12345Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(1+t*c12p + 3*(2+t*c345p)) # B phase
        
        fBGLRed = alphaRed*3*(delta**2) + 3*beta12345Red*(delta**4)

        beta245Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(2+t*c245p) # A Phase
        fAGLRed = alphaRed*2*(delta**2) + 4*beta245Red*(delta**4)

        fBGL_array[indexT,indexDelta] = fBGLRed
        fAGL_array[indexT,indexDelta] = fAGLRed

    print('fGLRed_Delta is:', fBGL_array[indexT,:])
    print('fGLRed_Delta is:', fAGL_array[indexT,:])
#    plot1.plot(Delta, fGL_array[indexT,:], 'o-'); plot1.ylabel('ev^2'); plot1.xlabel('ev')
#    plot1.show()
    

# Plot all data
for iT in range(0, lengthT, 1):
    #indexT = math.floor(T/stepT)
    indexT = iT
    plot1.plot(Delta, fBGL_array[indexT,:], '-.'); plot1.ylabel(r'$f^{B}_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$\Delta$/ev')
    plot1.plot(Delta, fAGL_array[indexT,:], '-'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$\Delta$/ev')

plot1.show()    
    



