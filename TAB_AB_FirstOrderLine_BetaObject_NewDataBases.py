# this script is for calculating the TAB, which is the equlibrium A-B 1st order line
# the strong coupling corrected \beta coefficients are taken from PRB. 92. 144515, 
# the priciple of calculation is solving equation \beta_{A} = \beta_{B}.
# the out put plot coincide with result in Fig.1 of PRB. 92. 144515 

# from math import pi

# zeta3 = 1.2020569;

# print(pi," ",zeta3);

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

pressure = np.arange(0.0, 34.0, 0.05);
print(" pressure arange is ", pressure)

C1=[];C2=[];C3=[];C4=[];C5=[];Tc=[]

TAB = [];

BetaObject = BETA('CiCoefficients')

for p in pressure:
      #p = pressure[i]
      print(" Now pressure is ",p)

      judgementList = region_judgement(p)

      print('judgementList is', judgementList)
      print('low pressure is: ',judgementList[0], ',high pressure is: ',judgementList[1], ',interpolation region is: ', judgementList[2])

      pk = judgementList[0]; pk1 = judgementList[1]; k = judgementList[2]


      

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

      tab = (tcp/(c1p + c2p + (1/3)*(c3p + c4p + c5p)-(c2p + c4p + c5p)))*(1/3)
      TAB.append(tab)

      print('\n\n tab now is: ',tab,'\n\n\n' );
      

     
#     print(c1[i], c3[i], c4[i], c5[i], Tc[i]);


#     
#     Tab = tab * Tc[i];
#     TAB.append(Tab); 
#     # print(Tab);

#for i in range(0, 18):
#    print(' for pressure ', p[i],'bar',' Tab is ',TAB[i]);

plot1.plot(pressure, TAB, 'o-');plot1.xlabel('p/bar'); plot1.ylabel('T_{AB}/mK')
plot1.savefig('TAB_p.pdf')
plot1.clf()
plot1.cla()
plot1.close()

plot1.plot(TAB, pressure, 'o-'); plot1.ylabel('p/bar'); plot1.xlabel('T_{AB}/mK');
plot1.savefig('p_TAB.pdf');
plot1.clf()
plot1.cla()
plot1.close()


startpoint = round(21.0/0.05);endpoint = round(34.0/0.05);
plot1.plot(TAB[startpoint:endpoint], pressure[startpoint:endpoint], 'o-'); plot1.ylabel('p/bar'); plot1.xlabel('T_{AB}/mK');
plot1.savefig('p_TAB.pdf');
#plot1.show()
plot1.clf()
plot1.cla()
plot1.close()

plot1.plot(Tc, pressure, 'o-');plot1.xlabel('p/bar'); plot1.ylabel('T_{AB}/mK')
startpoint = round(21.1/0.05);endpoint = round(34.0/0.05);
plot1.plot(TAB[startpoint:endpoint], pressure[startpoint:endpoint], 'o-'); plot1.ylabel('p/bar'); plot1.xlabel('T_{AB}/mK');
plot1.savefig('p_TAB_TC.pdf');
plot1.show()




