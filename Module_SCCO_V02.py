################################################################################################
#############              Introduction and Description                            #############
################################################################################################
#'''
#   >>>> READ ME! Please!
#'''
# This script is the class module of intepolations of strong coupling (SC) correction data,
# such as Tc, effective mass, Fermi velocity. They are necessary for calculating the SC corrected 
# GL coefficients e.g., \alpha and \beta_{i}. The later is also called " material parameters ".

# caution!
#'''
#  Users who focus on the SC corrected paramters e.g., gap energies of A/B phase, material paramters
#  are NOT recommended to use this class module directly. 
  
#  There is other module named " Module_SC_Beta_V*.py " is recommened to use in this case.
  
#'''

# The data sheet comes form table. II of
# PRB 101.024517

# This is version 0.1 (v.0.1) of the intepolation object (class)

# author: quang. zhang (github@timohyva)

# zeta3 = 1.2020569
  
#################################################################################################  

import numpy as np

# data sheet from talbe.II of PRB. 101. 024517
P = np.arange(0.0, 36.0, 2.0) # pressure date, step 2.0 bar
c1 = [-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413]
c2 = [-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645]
c3 = [-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268]
c4 = [-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518]
c5 = [-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815]
Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]
Ms = [2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82] # in unit of helium-3 atom
VF = [59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23]
XI0 = [77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76] # in unit of nm (10^-9 m)
         
class BETA:

     # class of beta object
    

     def __init__(self,name):
         self.name = name
         print(" beta oject is crated ")

     def c1_function(self,P,c1,pressure):
        self.c1p =  np.interp(pressure,P,c1)
    
     def c2_function(self,P,c2,pressure):
         self.c2p =  np.interp(pressure,P,c2)

     def c3_function(self,P,c3,pressure):
         self.c3p =  np.interp(pressure,P,c3)

     def c4_function(self,P,c4,pressure):
         self.c4p =  np.interp(pressure,P,c4)

     def c5_function(self,P,c5,pressure):
         self.c5p =  np.interp(pressure,P,c5)

     def tc_function(self,P,Tc,pressure):
         self.tcp =  np.interp(pressure,P,Tc)

     def mStar_function(self,P,MS,pressure):
         self.ms =  np.interp(pressure,P,MS)

     def vFermi_function(self,P,VF,pressure):
         self.vf = np.interp(pressure,P,VF)
         
     def xi0_function(self,P,XI0,pressure):
         self.xi0 = np.interp(pressure,P,XI0)    
         
         



