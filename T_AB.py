# this script is for calculating the TAB, which is the equlibrium A-B 1st order line
# the strong coupling corrected \beta coefficients are taken from PRB. 92. 144515, 
# the priciple of calculation is solving equation \beta_{A} = \beta_{B}.
# the out put plot coincide with result in Fig.1 of PRB. 92. 144515 

# from math import pi

# zeta3 = 1.2020569;

# print(pi," ",zeta3);

import matplotlib.pyplot as ploot



c1 = [0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, -0.01, -0.01, -0.01, -0.02, -0.02, -0.02, -0.03, -0.03];
c3 = [0.10, -0.14, -0.24, -0.28, -0.30, -0.31, -0.31, -0.30, -0.27, -0.27, -0.26, -0.26, -0.26, -0.27, -0.27, -0.28, -0.27, -0.27];
c4 = [-0.15, -0.37, -0.48, -0.54, -0.58, -0.60, -0.61, -0.62, -0.66, -0.68, -0.69, -0.71, -0.72, -0.73, -0.74, -0.74, -0.75, -0.75];
c5 = [0.16, 0.19, 0.19, 0.18, 0.17, 0.15, 0.13, 0.11, 0.10, 0.09, 0.07, 0.06, 0.04, 0.03, 0.01, -0.01, -0.02, -0.03];

Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486];

p = range(0, 36, 2);

TAB = [];

for i in range(0, 18):
     print(c1[i], c3[i], c4[i], c5[i], Tc[i]);


     tab = (1/3)/(c1[i] + (1/3)*c3[i] + (-2/3)*c4[i] + (-2/3)*c5[i]);
     Tab = tab * Tc[i];
     TAB.append(Tab); 
     # print(Tab);

for i in range(0, 18):
     print(' for pressure ', p[i],'bar',' Tab is ',TAB[i]);


ploot.plot(p[10:18], TAB[10:18], 'o-'); ploot.xlabel('p/bar'); ploot.ylabel('T_{AB}/mK');
ploot.show();



