# this script is for plot-check the GL free energy density under given pressure
# the strong coupling corrected \beta coefficients are taken from PRB. 92. 144515, 

# the out put plot coincide with result in Fig.1 of PRB. 92. 144515 



# zeta3 = 1.2020569;



import matplotlib.pyplot as plot1
import numpy as np

from math import pi




def intepolation_presure_SC(pk1,pk,p,fk1,fk):
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
         c1 = [0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, -0.01, -0.01, -0.01, -0.02, -0.02, -0.02, -0.03, -0.03]
         self.c1p =  intepolation_presure_SC(pk1,pk,pressure,c1[k+1],c1[k])

     def c2_function(self,pressure,pk,pk1,k):
         c2 = [-0.11, -0.04, -0.01, -0.01, -0.02, -0.03, -0.04, -0.05, -0.05, -0.06, -0.06, -0.07, -0.07, -0.07, -0.07, -0.07, -0.07, -0.07]
         self.c2p =  intepolation_presure_SC(pk1,pk,pressure,c2[k+1],c2[k])

     def c3_function(self,pressure,pk,pk1,k):
         c3 = [0.10, -0.14, -0.24, -0.28, -0.30, -0.31, -0.31, -0.30, -0.27, -0.27, -0.26, -0.26, -0.26, -0.27, -0.27, -0.28, -0.27, -0.27]
         self.c3p =  intepolation_presure_SC(pk1,pk,pressure,c3[k+1],c3[k])

     def c4_function(self,pressure,pk,pk1,k):
         c4 = [-0.15, -0.37, -0.48, -0.54, -0.58, -0.60, -0.61, -0.62, -0.66, -0.68, -0.69, -0.71, -0.72, -0.73, -0.74, -0.74, -0.75, -0.75]
         self.c4p =  intepolation_presure_SC(pk1,pk,pressure,c4[k+1],c4[k])

     def c5_function(self,pressure,pk,pk1,k):
         c5 = [0.16, 0.19, 0.19, 0.18, 0.17, 0.15, 0.13, 0.11, 0.10, 0.09, 0.07, 0.06, 0.04, 0.03, 0.01, -0.01, -0.02, -0.03]
         self.c5p =  intepolation_presure_SC(pk1,pk,pressure,c5[k+1],c5[k])

     def tc_function(self,pressure,pk,pk1,k):
         Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]
         self.tcp =  intepolation_presure_SC(pk1,pk,pressure,Tc[k+1],Tc[k])
         
         


# main code starts from here

# pressure = np.arange(0.0, 34.05, 0.05);

zeta3 = 1.2020569;
kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1

pressure = 24.0; #bar Unit
print(" pressure arange is ", pressure)
p = pressure

C1=[];C2=[];C3=[];C4=[];C5=[];Tc=[];

TAB = [];fGLRed_delta=[]

BetaObject = BETA('CiCoefficients')
#p = pressure[i]
print(" Now pressure is ",p)

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

print("pk = ",pk, "pk1 =", pk1, "k = ", k)   
      
BetaObject.c1_function(p,pk,pk1,k);c1p = BetaObject.c1p
print(" c1p =",c1p, " object.c1 = ",BetaObject.c1p)
      
BetaObject.c2_function(p,pk,pk1,k);c2p = BetaObject.c2p
BetaObject.c3_function(p,pk,pk1,k);c3p = BetaObject.c3p
BetaObject.c4_function(p,pk,pk1,k);c4p = BetaObject.c4p
BetaObject.c5_function(p,pk,pk1,k);c5p = BetaObject.c5p
BetaObject.tc_function(p,pk,pk1,k);tcp = BetaObject.tcp

print("c1p =",c1p,"c2p =",c2p,"c3p =", c3p,"c4p =", c4p, "c5p =",c5p,"tcp =", tcp)
      
C1.append(c1p)
C2.append(c2p)
C3.append(c3p)
C4.append(c4p)
C5.append(c5p)
Tc.append(tcp)



T = 0.0*(10**-3) #Kelvin
Delta = np.arange(0.0, 10000*kb*tcp*(10**(-3)), 10**(-6))

for delta in Delta:
      
      
      t = T/tcp
      alphaRed = (1/3)*(t-1);beta245Red = (1/(pi*pi*kb*kb))*(1/30)*(7/8)*zeta3*(2/(tcp*tcp)+(t/(tcp*tcp))*(c2p+c4p+c5p))
      fGLRed = alphaRed*2*(delta**2) + 4*beta245Red*(delta**4)

      fGLRed_delta.append(fGLRed)

      print("fGLRed",fGLRed)
      
# Plot the data
      
#plot1.plot(pressure, TAB, 'o-');plot1.xlabel('p/bar'); plot1.ylabel('T_{AB}/mK')
#plot1.savefig('TAB_p.pdf')
#plot1.clf()
#plot1.cla()
#plot1.close()

#startpoint = round(21.0/0.05);endpoint = round(34.0/0.05);
plot1.plot(Delta, fGLRed_delta, 'o-'); plot1.ylabel('ev^2'); plot1.xlabel('ev');
plot1.savefig('fGL_A.pdf');
#plot1.show()
plot1.clf()
plot1.cla()
plot1.close()





