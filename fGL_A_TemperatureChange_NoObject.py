# this script is for plot-check the GL free energy density under given pressure
# the strong coupling corrected \beta coefficients are taken from PRB. 92. 144515, 

# the out put plot coincide with result in Fig.1 of PRB. 92. 144515 



# zeta3 = 1.2020569;



import matplotlib.pyplot as plot1
import numpy as np

from math import pi
import math


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

pressure = 32.0; #bar Unit
print(" pressure arange is ", pressure)
p = pressure

C1=[];C2=[];C3=[];C4=[];C5=[];Tc=[];

TAB = [];fGLRed_delta=[]

print(" Now pressure is ",p)

c2p = -0.1583; c4p = -0.3388; c5p = -0.3717
c245p = c2p + c4p + c5p
Tcp = 2.463*(10**(-3))

print("c2p =",c2p,"c4p =", c4p, "c5p =",c5p,"tcp =", Tcp)
      
stepT = 0.25*(10**-3) 
Temperature = np.arange(0.0, 2.5*(10**-3), stepT) #Kelvin
#Temperature = [0.0]

stepDelta = 0.1*kb*Tcp
Delta = np.arange(0.0, 6*kb*Tcp, stepDelta)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Delta is',Delta,'\n length of Delta is ', len(Delta))
lengthDelta = len(Delta)

fGL_array = np.ones((lengthT,lengthDelta))

for T in Temperature:
    print('\n\n Now T is:', T, '\n\n')
    indexT = math.floor(T/stepT)
    print('indexT is ',indexT)
    
    for delta in Delta:
        indexDelta = math.floor(delta/stepDelta)
        print('indexDelta is',indexDelta)
       
        t = T/Tcp
        #print('pressure is:, ',p, 'Tcp is:,', tcp)
          
        alphaRed = (1/3)*(t-1)
        beta245Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(2+t*c245p)
        fGLRed = alphaRed*2*(delta**2) + 4*beta245Red*(delta**4)

        fGL_array[indexT,indexDelta] = fGLRed

        #print("fGLRed",fGLRed)
      

#startpoint = round(21.0/0.05);endpoint = round(34.0/0.05);
    print('fGLRed_Delta is:', fGL_array[indexT,:])
#    plot1.plot(Delta, fGL_array[indexT,:], 'o-'); plot1.ylabel('ev^2'); plot1.xlabel('ev')
#    plot1.show()
    fGLRed_delta=[] # clear old free energy

## Plot all data
for T in Temperature:
    indexT = math.floor(T/stepT)
    plot1.plot(Delta, fGL_array[indexT,:], '-'); plot1.ylabel('f_{GL}.N(o)⁻1/ev^2'); plot1.xlabel('\Delta/ev')

plot1.show()    
    



