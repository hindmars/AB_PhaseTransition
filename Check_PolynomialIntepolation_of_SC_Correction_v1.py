# This script compare the polynoimal intepolations of the strong coupling corrections
# of beta in the table. II and the numbers in the data sheet table.I of PRB. 92. 144515. 

# The reason to do this is the discrapancy between table.I and table. II in PRB. 101. 024517

# author: Quang. Zhang/timohyva:Github

#################################################
#################################################

# zeta3 = 1.2020569;



import matplotlib.pyplot as plot1
import numpy as np

from math import pi



def PolyNomial_intepolation_presure_SC(p,interpolation_Coes):
    a0 = interpolation_Coes[0]
    a1 = interpolation_Coes[1]
    a2 = interpolation_Coes[2]
    a3 = interpolation_Coes[3]
    a4 = interpolation_Coes[4]
    a5 = interpolation_Coes[5]
    a6 = interpolation_Coes[6]
    a7 = interpolation_Coes[7]
    a8 = interpolation_Coes[8]

    # polynimial interpolation in table.II of PRB. 92. 144515
    fip = a0 + a1*p + a2*p*p + a3*p*p*p + a4*p*p*p*p + a5*p*p*p*p*p + a6*p*p*p*p*p*p + a7*p*p*p*p*p*p*p + a8*p*p*p*p*p*p*p*p

    # return the cip or Tcp
    return fip
         
class BETA:

     # class of beta object
    

     def __init__(self,name):
         self.name = name
         print(" beta oject is crated ")

     def c1_function(self,pressure):
         
         interpolation_Coefficients = [3.070*10**(-2), -2.081*10**(-3), 2.133*10**(-5), -4.189*10**(-7), 0.0, 0.0, 0.0, 0.0, 0.0]
         self.c1p =  PolyNomial_intepolation_presure_SC(pressure,interpolation_Coefficients)

     def c2_function(self,pressure):

         interpolation_Coefficients = [-1.074*10**(-1), 5.412*10**(-2), -1.081*10**(-2), 1.025*10**(-3), -5.526*10**(-5), 1.722*10**(-6), -2.876*10**(-8), 1.991*10**(-10), 0.0]
         self.c2p =  PolyNomial_intepolation_presure_SC(pressure,interpolation_Coefficients)

     def c3_function(self,pressure):
         
         interpolation_Coefficients = [1.038*10**(-1), -1.752*10**(-1), 3.488*10**(-2), -4.243*10**(-3), 3.316*10**(-4), -1.623*10**(-5), 4.755*10**(-7), -7.587*10**(-9), 5.063*10**(-11)]
         self.c3p =  PolyNomial_intepolation_presure_SC(pressure,interpolation_Coefficients)

     def c4_function(self,pressure):
         
         interpolation_Coefficients = [-1.593*10**(-1), -1.350*10**(-1), 1.815*10**(-2), -1.339*10**(-3), 5.316*10**(-5), -1.073*10**(-6), 8.636*10**(-9), 0.0, 0.0]
         self.c4p =  PolyNomial_intepolation_presure_SC(pressure,interpolation_Coefficients)

     def c5_function(self,pressure):
         
         interpolation_Coefficients = [1.610*10**(-1), 2.263*10**(-2), -4.921*10**(-3), 3.810*10**(-4), -1.529*10**(-5), 3.071*10**(-7), -2.438*10**(-9), 0.0, 0.0]
         self.c5p =  PolyNomial_intepolation_presure_SC(pressure,interpolation_Coefficients)

     # def tc_function(self,pressure):
     #     Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]
     #     self.tcp =  intepolation_presure_SC(pk,pk1,pressure,Tc[k],Tc[k+1])
         
         
# the data sheet in table. I of PRB. 92. 144515
PressureSheet = np.arange(0.0, 36.0, 2.0)
c1_sheet = [0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, -0.01, -0.01, -0.01 , -0.02, -0.02, -0.02, -0.03, -0.03]
c2_sheet = [-0.11, -0.04, -0.01, -0.01, -0.02, -0.03, -0.04, -0.05, -0.05, -0.06, -0.06, -0.07, -0.07, -0.07, -0.07, -0.07, -0.07, -0.07]
c3_sheet = [0.10, -0.14, -0.24, -0.28, -0.30, -0.31, -0.31, -0.30, -0.27, -0.27, -0.26, -0.26, -0.26, -0.27, -0.27, -0.28, -0.27, -0.27]
c4_sheet = [-0.15, -0.37, -0.48, -0.54, -0.58, -0.60, -0.61, -0.62, -0.66, -0.68, -0.69, -0.71, -0.72, -0.73, -0.74, -0.74, -0.75, -0.75]
c5_sheet = [0.16, 0.19, 0.19, 0.18, 0.17, 0.15, 0.13, 0.11, 0.10, 0.09, 0.07, 0.06, 0.04, 0.03, 0.01, -0.01, -0.02, -0.03]
Cn_sheet = np.vstack((c1_sheet, c2_sheet, c3_sheet, c4_sheet, c5_sheet))
     
# test the object

BetaObject = BETA('betaAndTc')

# p = 32 # pressure, bar
pressureArray = np.arange(0.0, 34.1, 0.1)

# saving lists
CArray = np.zeros((5,len(pressureArray)))

# for p in pressureArray:
for Indp in range(0, len(pressureArray), 1):

    p = pressureArray[Indp]
    
    # judgementList = region_judgement(p)

    # print('judgementList is', judgementList)

    # print('low pressure is: ',judgementList[0], ',high pressure is: ',judgementList[1], ',interpolation region is: ', judgementList[2])

    # pk = judgementList[0]; pk1 = judgementList[1]; k = judgementList[2]

    BetaObject.c1_function(p);c1p = BetaObject.c1p;CArray[0,Indp] = c1p # C1.append(c1p)
    BetaObject.c2_function(p);c2p = BetaObject.c2p;CArray[1,Indp] = c2p # C2.append(c2p)
    BetaObject.c3_function(p);c3p = BetaObject.c3p;CArray[2,Indp] = c3p # C3.append(c3p)
    BetaObject.c4_function(p);c4p = BetaObject.c4p;CArray[3,Indp] = c4p # C4.append(c4p)
    BetaObject.c5_function(p);c5p = BetaObject.c5p;CArray[4,Indp] = c5p # C5.append(c5p)
    # BetaObject.tc_function(p);Tcp = BetaObject.tcp


    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p is ',c5p,'\n\n\n')


# print(CArray, '\n\n')
print(Cn_sheet)


# Plot polynomial intepolation data and sheet in table.II of PRB. 92. 144515,
# The sequence of plot-showing is from c1 till c5 (beta1 till beta5), dots plots are table.I in sheet.
for ii in range(0, 5, 1):
    plot1.plot(pressureArray, CArray[ii,:], '-'); # plot1.ylabel(r'$f^{B}_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$p$/bar')
    # plot1.plot(Delta, fAGL_array[indexT,:], '-'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$p$/bar')
    plot1.plot(PressureSheet, Cn_sheet[ii,:], 'o'); plot1.ylabel(r'SC Coefficient'); plot1.xlabel(r'$p$/bar')
    plot1.show()     
