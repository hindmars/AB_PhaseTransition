#############################################################
#############################################################
'''Module for generating p, T arries from different experimental paths.

There are three types of path: Cons-p, Cons-Q and other according to Cornell's data

author: Quang. Zhang (timohyva@github)
'''


import numpy as np
import pandas as pd

import Module_SC_Beta_V05 as SC_beta
import he3_tools_Vn01 as hn


#############################################################
##                    turn on PLTS                         ##
#############################################################

SC_beta.turn_on_PLTS()

#############################################################
##                Constant-p data source                   ##
#############################################################

# import the csv data of Parpia's new data
sheet_consp = np.genfromtxt('Cornell22_const_P.csv', skip_header=1, delimiter=',')
# pressure

# pressure
pressure = sheet_consp[:, 0];
# print("\n\n events' pressures in constant-p runs:\n ",pressure)

p_constp_arr = np.zeros(pressure.shape)

# T_IC, unit mK
# T_IC = sheet_consp[:, 2]; print("\n T_IC events' temperatures in constant-p runs:\n ", T_IC)

# T_HEC
T_HEC = sheet_consp[:, 2];
# print("\n T_HEC events' temperatures in constant-p runs:\n ", T_HEC)

# constrct the TAB array of constant-p
TAB_constp = np.zeros(pressure.shape); # print("\n TAB_constp looks like ", TAB_constp)
Tc_constp = np.zeros(pressure.shape)

# print(pressure[(pressure != np.nan)])
for ip in np.arange(0,len(pressure), 1):
  if np.isnan(SC_beta.TAB_RWSco(pressure[ip])):
   continue
  else:
    # print("\n ", SC_beta.TAB_RWSco(pressure[ip]))
    TAB_constp[ip] = SC_beta.TAB_RWSco(pressure[ip])
    Tc_constp[ip] = SC_beta.Tcp(pressure[ip])
    p_constp_arr[ip] = pressure[ip]
    
# print(" \nTAB_constp looks like\n ",TAB_constp, " \nafter filting\n ",TAB_constp[TAB_constp != 0.])
# print(" \nTc_constp looks like\n ",Tc_constp, " \nafter filting\n ",Tc_constp[Tc_constp != 0.])
# print(" \np_arr looks like\n ", p_constp_arr, " \nafter filting\n ",p_constp_arr[p_constp_arr != 0.])

# convert Temperatures to PLTS scate
# TAB_constp = hn.T_GtoPlts6_low_poly(1000*TAB_constp[TAB_constp != 0.])
# Tc_constp = hn.T_GtoPlts6_low_poly(1000*Tc_constp[Tc_constp != 0.])
# p_constp_arr = p_constp_arr[p_constp_arr != 0.]; # print(" \n length of p_constp_arr ", len(p_constp_arr))


###############################################################
##                 constant Q data source                    ##
###############################################################


##   >>>>>>>>>>>>>>>>   Q1 data     <<<<<<<<<<<<<<<<<        ##

stp1 = 0.001
T1arr1 = np.arange(2.058, 2.25, stp1)
T1arr2 = np.arange(2.25, stp1+SC_beta.Tcp(29.3)*1000, stp1)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ1 = 25.370095287443853; TABQ1 = 2.143458570116102 # mK

T1_HEC_arr = np.append(T1arr1, T1arr2)
boolen_Q1 = T1_HEC_arr > TABQ1
T1_HEC_arr = T1_HEC_arr[boolen_Q1]

p1_arr = np.append(np.linspace(22.0626, 29.3, len(T1arr1)), 29.3*np.ones(T1arr2.shape))
p1_arr = p1_arr[boolen_Q1]

# print(" \nlength of p1_arr: ",len(p1_arr), " length of T1_HEC_arr: ", len(T1_HEC_arr))                                                      

##   >>>>>>>>>>>>>>>>   Q2 data     <<<<<<<<<<<<<<<<<        ##

stp2 = 0.001
T2arr1 = np.arange(2.102, 2.275, stp2)
T2arr2 = np.arange(2.275, stp2+SC_beta.Tcp(29.3)*1000, stp1)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ2 = 24.43600969915881; TABQ2 = 2.1677548560249678 # mK
                                                          
T2_HEC_arr = np.append(T2arr1, T2arr2)
boolen_Q2 = T2_HEC_arr > TABQ2
T2_HEC_arr = T2_HEC_arr[boolen_Q2]

p2_arr = np.append(np.linspace(21.5357, 29.3, len(T2arr1)), 29.3*np.ones(T2arr2.shape))
p2_arr = p2_arr[boolen_Q2]

# print(" \nlength of p2_arr: ",len(p2_arr), " length of T2_HEC_arr: ", len(T2_HEC_arr))

##   >>>>>>>>>>>>>>>>   Q3 data     <<<<<<<<<<<<<<<<<        ##

stp3 = 0.001
T3arr1 = np.arange(2.084, 2.285, stp3)
T3arr2 = np.arange(2.285, stp3+SC_beta.Tcp(29.3)*1000, stp3)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ3 =  24.64062849593238; TABQ3 = 2.162728479519771 # mK
                                                          
T3_HEC_arr = np.append(T3arr1, T3arr2)
boolen_Q3 = T3_HEC_arr > TABQ3
T3_HEC_arr = T3_HEC_arr[boolen_Q3]

p3_arr = np.append(np.linspace(21.696, 29.3, len(T3arr1)), 29.3*np.ones(T3arr2.shape))
p3_arr = p3_arr[boolen_Q3]

# print(" \nlength of p3_arr: ",len(p3_arr), " length of T3_HEC_arr: ", len(T3_HEC_arr))

##   >>>>>>>>>>>>>>>>   Q4 data     <<<<<<<<<<<<<<<<<        ##

stp4 = 0.001
T4arr1 = np.arange(2.121, 2.29, stp4)
T4arr2 = np.arange(2.29, stp4+SC_beta.Tcp(29.3)*1000, stp4)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ4 =  23.857740367654966; TABQ4 = 2.1819963137662786 # mK
                                                          
T4_HEC_arr = np.append(T4arr1, T4arr2)
boolen_Q4 = T4_HEC_arr > TABQ4
T4_HEC_arr = T4_HEC_arr[boolen_Q4]

p4_arr = np.append(np.linspace(20.919, 29.3, len(T4arr1)), 29.3*np.ones(T4arr2.shape))
p4_arr = p4_arr[boolen_Q4]

# print(" \nlength of p4_arr: ",len(p4_arr), " length of T4_HEC_arr: ", len(T4_HEC_arr))


##    >>>>>>>>>>>>>>>>   Q5 data, Q 65    <<<<<<<<<<<<<<<<<<     ##

stp5 = 0.001
T5arr1 = np.arange(2.146, 2.306, stp5)
T5arr2 = np.arange(2.306, stp5+SC_beta.Tcp(29.3)*1000, stp5)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ5 = 23.267794143331198; TABQ5= 2.196113989637306 # mK
                                                          
T5_HEC_arr = np.append(T5arr1, T5arr2)
boolen_Q5 = T5_HEC_arr > TABQ5
T5_HEC_arr = T5_HEC_arr[boolen_Q5]

p5_arr = np.append(np.linspace(20.553, 29.3, len(T5arr1)), 29.3*np.ones(T5arr2.shape))
p5_arr = p5_arr[boolen_Q5]

# print(" \nlength of p5_arr: ",len(p5_arr), " length of T5_HEC_arr: ", len(T5_HEC_arr))

##    >>>>>>>>>>>>>>>>   Q6 data,Q 53    <<<<<<<<<<<<<<<<<<     ##

stp6 = 0.001
T6arr1 = np.arange(2.181, 2.323, stp6)
T6arr2 = np.arange(2.323, stp6+SC_beta.Tcp(29.3)*1000, stp6)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ6 = 22.46725994061827; TABQ6= 2.212642487046632 # mK
                                                          
T6_HEC_arr = np.append(T6arr1, T6arr2)
boolen_Q6 = T6_HEC_arr > TABQ6
T6_HEC_arr = T6_HEC_arr[boolen_Q6]

p6_arr = np.append(np.linspace(20.47, 29.3, len(T6arr1)), 29.3*np.ones(T6arr2.shape))
p6_arr = p6_arr[boolen_Q6]

# print(" \nlength of p6_arr: ",len(p6_arr), " length of T6_HEC_arr: ", len(T6_HEC_arr))


##    >>>>>>>>>>>>>>>>   Q7 data, Q 45    <<<<<<<<<<<<<<<<<<     ##

stp7 = 0.001
T7arr1 = np.arange(2.184, 2.339, stp7)
T7arr2 = np.arange(2.339, stp7+SC_beta.Tcp(29.3)*1000, stp7)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ7 = 22.127622984223088; TABQ7 = 2.218911917098446 # mK

                                                          
T7_HEC_arr = np.append(T7arr1, T7arr2)
boolen_Q7 = T7_HEC_arr > TABQ7
T7_HEC_arr = T7_HEC_arr[boolen_Q7]

p7_arr = np.append(np.linspace(20.025, 29.3, len(T7arr1)), 29.3*np.ones(T7arr2.shape))
p7_arr = p7_arr[boolen_Q7]

# print(" \nlength of p7_arr: ",len(p7_arr), " length of T7_HEC_arr: ", len(T7_HEC_arr))


##    >>>>>>>>>>>>>>>>   Q8 data, Q 40    <<<<<<<<<<<<<<<<<<     ##

stp8 = 0.001
T8arr1 = np.arange(2.198, 2.351, stp8)
T8arr2 = np.arange(2.351, stp8+SC_beta.Tcp(29.3)*1000, stp8)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ8 = 21.836500203760842; TABQ8 = 2.2240414507772024 # mK

                                                          
T8_HEC_arr = np.append(T8arr1, T8arr2)
boolen_Q8 = T8_HEC_arr > TABQ8
T8_HEC_arr = T8_HEC_arr[boolen_Q8]

p8_arr = np.append(np.linspace(20.255, 29.3, len(T8arr1)), 29.3*np.ones(T8arr2.shape))
p8_arr = p8_arr[boolen_Q8]

# print(" \nlength of p8_arr: ",len(p8_arr), " length of T8_HEC_arr: ", len(T8_HEC_arr))

##    >>>>>>>>>>>>>>>>   Q9 data, " const Q, waiting one day then cooling "  <<<<<<<<<<<<<<<<<<     ##

stp9 = 0.001
T9arr1 = np.arange(2.146, 2.306, stp9)
T9arr2 = np.arange(2.306, stp9+SC_beta.Tcp(29.3)*1000, stp9)
# print(" \nT1arr1 looks like\n ", T1arr1)
# print(" \nT1arr2 looks like\n ", T1arr2)

# TAb point digtized from the plot in cornell's manuscript 
pABQ9 = 23.267794143331198; TABQ9 = 2.196113989637306 # mK
                                                          
T9_HEC_arr = np.append(T9arr1, T9arr2)
boolen_Q9 = T9_HEC_arr > TABQ9
T9_HEC_arr = T9_HEC_arr[boolen_Q9]

p9_arr = np.append(np.linspace(20.553, 29.3, len(T9arr1)), 29.3*np.ones(T9arr2.shape))
p9_arr = p9_arr[boolen_Q9]

# print(" \nlength of p9_arr: ",len(p9_arr), " length of T9_HEC_arr: ", len(T9_HEC_arr))

Tarr_constQ_list = [T1_HEC_arr, T2_HEC_arr, T3_HEC_arr, T4_HEC_arr, T5_HEC_arr, T6_HEC_arr, T7_HEC_arr, T8_HEC_arr, T9_HEC_arr]
parr_constQ_list = [p1_arr, p2_arr, p3_arr, p4_arr, p5_arr, p6_arr, p7_arr, p8_arr, p9_arr]


###############################################################
##               others source group no. 1                   ##
###############################################################

##   >>>>>>>>>>>>>>>>   O1_1 data     <<<<<<<<<<<<<<<<<      ##

stpO1 = 0.001
TO1arr1 = np.arange(2.012, 2.14, stpO1)
TO1arr2 = 2.14*np.ones(100)
TO1arr3 = np.arange(2.14+stpO1, stpO1+SC_beta.Tcp(29.3)*1000, stpO1)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO11 =  25.51887445887446; TABO11 = 2.1395579268292684 # mK
                                                         
TO1_HEC_arr = np.append(np.append(TO1arr1, TO1arr2), TO1arr3)
boolen_O11 = TO1_HEC_arr > TABO11
TO1_HEC_arr = TO1_HEC_arr[boolen_O11]


pO1_arr = np.append(np.append(23.*np.ones(TO1arr1.shape), np.linspace(23., 29.3, len(TO1arr2))), 29.3*np.ones(TO1arr3.shape))
pO1_arr = pO1_arr[boolen_O11]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))

##   >>>>>>>>>>>>>>>>   O1_2 data     <<<<<<<<<<<<<<<<<      ##

stpO2 = 0.001
TO2arr1 = np.arange(1.949, 2.05, stpO2)
TO2arr2 = 2.05*np.ones(100)
TO2arr3 = np.arange(2.05+stpO2, stpO2+SC_beta.Tcp(29.3)*1000, stpO2)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO12 =  28.782770562770562; TABO12 = 2.0536204268292684 # mK
                                                          
TO2_HEC_arr = np.append(np.append(TO2arr1, TO2arr2), TO2arr3)
boolen_O12 = TO2_HEC_arr > TABO12
TO2_HEC_arr = TO2_HEC_arr[boolen_O12]

pO2_arr = np.append(np.append(25.*np.ones(TO2arr1.shape), np.linspace(25., 29.3, len(TO2arr2))), 29.3*np.ones(TO2arr3.shape))
pO2_arr = pO2_arr[boolen_O12]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))             

##   >>>>>>>>>>>>>>>>   O1_3 data     <<<<<<<<<<<<<<<<<      ##

stpO3 = 0.001
TO3arr1 = np.arange(2.123, 2.15, stpO3)
TO3arr2 = 2.15*np.ones(100)
TO3arr3 = np.arange(2.15+stpO3, stpO3+SC_beta.Tcp(29.3)*1000, stpO3)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO13 =  25.055151515151515; TABO13 = 2.1517149390243904 # mK
                                                          
TO3_HEC_arr = np.append(np.append(TO3arr1, TO3arr2), TO3arr3)
boolen_O13 = TO3_HEC_arr > TABO13
TO3_HEC_arr = TO3_HEC_arr[boolen_O13]


pO3_arr = np.append(np.append(21.*np.ones(TO3arr1.shape), np.linspace(21., 29.3, len(TO3arr2))), 29.3*np.ones(TO3arr3.shape))
pO3_arr = pO3_arr[boolen_O13]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))             

##   >>>>>>>>>>>>>>>>   O1_4 data     <<<<<<<<<<<<<<<<<      ##

stpO4 = 0.001
# TO4arr1 = np.arange(1.949, 2.05, stpO2)
TO4arr2 = 1.951*np.ones(100)
TO4arr3 = np.arange(1.951+stpO4, stpO4+SC_beta.Tcp(29.3)*1000, stpO4)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO14 =  29.3; TABO14 = 2.039786585365854 # mK
                                                          
TO4_HEC_arr = np.append(TO4arr2, TO4arr3)
boolen_O14 = TO4_HEC_arr > TABO14
TO4_HEC_arr = TO4_HEC_arr[boolen_O14]

pO4_arr = np.append(np.linspace(24.865, 29.3, len(TO4arr2)), 29.3*np.ones(TO4arr3.shape))
pO4_arr = pO4_arr[boolen_O14]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))

##   >>>>>>>>>>>>>>>>   O1_5 data     <<<<<<<<<<<<<<<<<      ##

stpO5 = 0.001
# TO4arr1 = np.arange(1.949, 2.05, stpO2)
TO5arr2 = 1.901*np.ones(100)
TO5arr3 = np.arange(1.901+stpO5, stpO5+SC_beta.Tcp(29.3)*1000, stpO5)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO15 =  29.3; TABO15 = 2.039786585365854 # mK

TO5_HEC_arr = np.append(TO5arr2, TO5arr3)
boolen_O15 = TO5_HEC_arr > TABO15
TO5_HEC_arr = TO5_HEC_arr[boolen_O15]

pO5_arr = np.append(np.linspace(27., 29.3, len(TO5arr2)), 29.3*np.ones(TO5arr3.shape))
pO5_arr = pO5_arr[boolen_O15]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))

##   >>>>>>>>>>>>>>>>   O1_6 data     <<<<<<<<<<<<<<<<<      ##

stpO6 = 0.001
TO6arr1 = np.arange(1.903, 2., stpO6)
TO6arr2 = 2.*np.ones(100)
TO6arr3 = np.arange(2.+stpO6, stpO6+SC_beta.Tcp(29.3)*1000, stpO6)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO16 =  29.3; TABO16 = 2.039786585365854 # mK
                                                          
TO6_HEC_arr = np.append(np.append(TO6arr1, TO6arr2), TO6arr3)
boolen_O16 = TO6_HEC_arr > TABO16
TO6_HEC_arr = TO6_HEC_arr[boolen_O16]


pO6_arr = np.append(np.append(27.*np.ones(TO6arr1.shape), np.linspace(27., 29.3, len(TO6arr2))), 29.3*np.ones(TO6arr3.shape))
pO6_arr = pO6_arr[boolen_O16]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))             

##   >>>>>>>>>>>>>>>>   O1_7 data     <<<<<<<<<<<<<<<<<      ##

stpO7 = 0.001
TO7arr1 = np.arange(1.959, 2.08, stpO7)
TO7arr2 = 2.08*np.ones(100)
TO7arr3 = np.arange(2.08+stpO7, stpO7+SC_beta.Tcp(29.3)*1000, stpO7)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)

# TAb point digtized from the plot in cornell's manuscript 
pABO17 =  27.6591341991342; TABO17 = 2.0829649390243903 # mK
                                                          
TO7_HEC_arr = np.append(np.append(TO7arr1, TO7arr2), TO7arr3)
boolen_O17 = TO7_HEC_arr > TABO17
TO7_HEC_arr = TO7_HEC_arr[boolen_O17]

pO7_arr = np.append(np.append(24.*np.ones(TO7arr1.shape), np.linspace(24., 29.3, len(TO7arr2))), 29.3*np.ones(TO7arr3.shape))
pO7_arr = pO7_arr[boolen_O17]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))

Tarr_other1_list = [TO1_HEC_arr, TO2_HEC_arr, TO3_HEC_arr, TO4_HEC_arr, TO5_HEC_arr, TO6_HEC_arr, TO7_HEC_arr]
parr_other1_list = [pO1_arr, pO2_arr, pO3_arr, pO4_arr, pO5_arr, pO6_arr, pO7_arr]

###############################################################
##               others source group no. 2                   ##
###############################################################

#    >>>>>>>>>>>>>>>>>>>>   O2-1   <<<<<<<<<<<<<<<<<<<<<     ##

stpO21 = 0.0002
# TO21arr1 = np.arange(2.0106, SC_beta.TAB_RWSco(23.)*1000, stpO21)
TO21arr1 = np.arange(2.0106, 2.203, stpO21)
# TO21arr2 = (SC_beta.TAB_RWSco(23.)*1000)*np.ones(200)
TO21arr2 = (2.203)*np.ones(200)
# TO21arr3 = np.arange(SC_beta.TAB_RWSco(23.)*1000+stpO21, stpO21+SC_beta.Tcp(23.)*1000, stpO21)
TO21arr3 = np.arange(2.203+stpO21, stpO21+SC_beta.Tcp(23.)*1000, stpO21)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)
                                                          
TO21_HEC_arr = np.append(np.append(TO21arr1, TO21arr2), TO21arr3)
boolen_O21 = ~(TO21_HEC_arr < (SC_beta.TAB_RWSco(23.)*1000))
TO21_HEC_arr = TO21_HEC_arr[boolen_O21]

pO21_arr = np.append(np.append(23.*np.ones(TO21arr1.shape), np.linspace(23., 27.5, len(TO21arr2))), 23.*np.ones(TO21arr3.shape))
pO21_arr = pO21_arr[boolen_O21]

# print(" \nlength of p1_arr: ",len(pO1_arr), " length of T1_HEC_arr: ", len(TO1_HEC_arr))

#    >>>>>>>>>>>>>>>>>>>>   O2-2   <<<<<<<<<<<<<<<<<<<<<     ##

stpO22 = 0.0002
# TO22arr1 = np.arange(2.0617, SC_beta.TAB_RWSco(23.)*1000, stpO22)
TO22arr1 = np.arange(2.0617, 2.203, stpO22)
# TO22arr2 = (SC_beta.TAB_RWSco(23.)*1000)*np.ones(200)
TO22arr2 = (2.203)*np.ones(200)
# TO22arr3 = np.arange(SC_beta.TAB_RWSco(23.)*1000+stpO22, stpO22+SC_beta.Tcp(23.)*1000, stpO22)
TO22arr3 = np.arange(2.203+stpO22, stpO22+SC_beta.Tcp(23.)*1000, stpO22)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)
                                                          
TO22_HEC_arr = np.append(np.append(TO22arr1, TO22arr2), TO22arr3)
boolen_O22 = ~(TO22_HEC_arr < (SC_beta.TAB_RWSco(23.)*1000))
TO22_HEC_arr = TO22_HEC_arr[boolen_O22]

pO22_arr = np.append(np.append(23.*np.ones(TO22arr1.shape), np.linspace(23., 26., len(TO22arr2))), 23.*np.ones(TO22arr3.shape))
pO22_arr = pO22_arr[boolen_O22]


#    >>>>>>>>>>>>>>>>>>>>   O2-3   <<<<<<<<<<<<<<<<<<<<<     ##

stpO23 = 0.0002
# TO23arr1 = np.arange(2.0996, SC_beta.TAB_RWSco(23.)*1000, stpO23)
TO23arr1 = np.arange(2.0996, 2.203, stpO23)
# TO23arr2 = (SC_beta.TAB_RWSco(23.)*1000)*np.ones(200)
TO23arr2 = (2.203)*np.ones(200)
# TO23arr3 = np.arange(SC_beta.TAB_RWSco(23.)*1000+stpO23, stpO23+SC_beta.Tcp(23.)*1000, stpO23)
TO23arr3 = np.arange(2.203+stpO23, stpO23+SC_beta.Tcp(23.)*1000, stpO23)
# print(" \nTO1arr1 looks like\n ", TO1arr1)
# print(" \nTO1arr2 looks like\n ", TO1arr2)
# print(" \nTO1arr3 looks like\n ", TO1arr3)
                                                          
TO23_HEC_arr = np.append(np.append(TO23arr1, TO23arr2), TO23arr3)
boolen_O23 = ~(TO23_HEC_arr < (SC_beta.TAB_RWSco(23.)*1000))
TO23_HEC_arr = TO23_HEC_arr[boolen_O23]

pO23_arr = np.append(np.append(23.*np.ones(TO23arr1.shape), np.linspace(23., 24.5, len(TO23arr2))), 23.*np.ones(TO23arr3.shape))
pO23_arr = pO23_arr[boolen_O23]

Tarr_other2_list = [TO21_HEC_arr, TO22_HEC_arr, TO23_HEC_arr]
parr_other2_list = [pO21_arr, pO22_arr, pO23_arr]


###############################################################
##                         no name                          ###
###############################################################

##   >>>>>>>>>>>>>>>>   NN_1 data     <<<<<<<<<<<<<<<<<      ##

stpNN1 = 0.001
TNN1arr1 = np.arange(2.059, 2.16, stpNN1)
TNN1arr2 = np.arange(2.16, stpNN1+SC_beta.Tcp(29.3)*1000, stpNN1)
# print(" \nTNN1arr1 looks like\n ", TNN1arr1)
# print(" \nTNN1arr3 looks like\n ", TNN1arr3)

# TAb point digtized from the plot in cornell's manuscript 
# pABNN1 =  29.3; TABNN1 = 2.039786585365854 # mK
                                                          
TNN1_HEC_arr = np.append(TNN1arr1, TNN1arr2)
boolen_NN1 = ~(TNN1_HEC_arr < (SC_beta.TAB_RWSco(29.3)*1000))
TNN1_HEC_arr = TNN1_HEC_arr[boolen_NN1]

pNN1_arr = np.append(np.linspace(28.579, 29.3, len(TNN1arr1)), 29.3*np.ones(TNN1arr2.shape))
pNN1_arr = pNN1_arr[boolen_NN1]

##   >>>>>>>>>>>>>>>>   NN_2 data     <<<<<<<<<<<<<<<<<      ##

stpNN2 = 0.001
TNN2arr1 = np.arange(2.109, 2.203, stpNN2)
TNN2arr2 = np.linspace(2.203, 2.203, 100)
TNN2arr3 = np.arange(2.203, stpNN2+SC_beta.Tcp(23)*1000, stpNN2)
                                                          
TNN2_HEC_arr = np.concatenate((TNN2arr1, TNN2arr2, TNN2arr3))
boolen_NN2 = ~(TNN2_HEC_arr < (SC_beta.TAB_RWSco(25)*1000))
TNN2_HEC_arr = TNN2_HEC_arr[boolen_NN2]

pNN2_arr = np.concatenate((np.linspace(25, 25, len(TNN2arr1)), 
                     np.linspace(25, 23, len(TNN2arr2)),
                     23*np.ones(TNN2arr3.shape)))
pNN2_arr = pNN2_arr[boolen_NN2]

Tarr_NN_list = [TNN1_HEC_arr, TNN2_HEC_arr]
parr_NN_list = [pNN1_arr, pNN2_arr]


###############################################################
##                      interfaces                           ##
###############################################################

def arr_generator_noname(no, Tarrl = Tarr_NN_list, parrl = parr_NN_list):
  '''Interface for handling array of noname run.

  *no* could be 1, 2
  '''
  # print(" number of data groups: ", no,", and len of Tarrl ", len(Tarrl))
  return parrl[no-1], Tarrl[no-1]


def arr_generator_other_2(no, Tarrl = Tarr_other2_list, parrl = parr_other2_list):
  '''Interface for handling arries generations for Other-1 

  *no* represnets the series number of datas. no could be 1,2,3.
  '''
  # print(" \n len of Tarrl ", len(Tarrl))
  return parrl[no-1], Tarrl[no-1]


def arr_generator_other_1(no, Tarrl = Tarr_other1_list, parrl = parr_other1_list):
  '''Interface for handling arries generations for Other-1 

  *no* represnets the series number of datas. no could be 1,2,3,4,5,6,7.
  '''
  # print(" \n len of Tarrl ", len(Tarrl))
  return parrl[no-1], Tarrl[no-1]

def arr_generator_constQ(no, Tarrl = Tarr_constQ_list, parrl = parr_constQ_list):
  '''Interface for handling arries generations for Const-Q

  *no* represnets the series number of datas. no could be 1,2,3,4,5,6,7,8.
  '''
  return parrl[no-1], Tarrl[no-1]
  

def arr_generator_constp(no, TAB_arr = TAB_constp, Tc_arr = Tc_constp, p_Arr = p_constp_arr):
  '''Interface for handling arries generations for Const-p

  *no* represnets the series number of datas. no couble be one integal between 1 and 39
  '''
  length = 200; 
  T_arr = 1000*np.linspace(TAB_arr[no-1], Tc_arr[no-1], length)
  p_arr = p_Arr[no-1]*np.ones(T_arr.shape)

  return p_arr, T_arr


def get_data(key, no):
  '''Interface for p, T arries generations

  calling: RR_pT.get_data(key, no), where key is one of following strings, 

      "constp", "constQ", "others_1(2)", "noname" and no is a int type number.
  
  return: p_arr, T_arr, T_arr in unit of mK
  
  For using this interface, module named " Module_SC_Beta_V05 " is necessary.
  '''
  if key == "constp": return arr_generator_constp(no)

  if key == "constQ": return arr_generator_constQ(no)

  if key == "others_1": return arr_generator_other_1(no)

  if key == "others_2": return arr_generator_other_2(no)

  if key == "noname": return arr_generator_noname(no)


def get_data_nos(key):
    """Counts number of paths in the dataset given by key.
    """
    if key == "others_2":
        return len(parr_other2_list)
    if key == "others_1":
        return len(parr_other1_list)
    if key == "constQ":
        return len(parr_constQ_list)
    if key == "constp":
        return len(p_constp_arr)
    if key == "noname":
        return len(parr_NN_list)
    return np.nan

