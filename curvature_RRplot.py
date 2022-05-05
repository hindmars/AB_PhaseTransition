#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 19:19:48 2022

plotting the curvature 1/R_crit of Cornell's experiments

const-p, const-Q, other-1 and other-2 are evaluated.

@author: Quang (timohyva@github)
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h

import Module_RR_pT as rr_pT
import Module_SC_Beta_V05 as SC_beta

import Module_New_Parpia_pT_parameter as Parpia

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")

################################################################
##                       cosnt-p data                         ##
################################################################


p_arr_constp = Parpia.pressure
T_IC_constp = Parpia.T_IC; t_IC_constp = T_IC_constp/(1000*SC_beta.Tcp(p_arr_constp))
T_HEC_constp =Parpia.T_HEC; t_HEC_constp = T_HEC_constp/(1000*SC_beta.Tcp(p_arr_constp))

Rcrit_HEC_constp = h.critical_radius(t_HEC_constp, p_arr_constp)
Rcrit_IC_constp = h.critical_radius(t_IC_constp, p_arr_constp)

print(" \n Rcrit_HEC ", Rcrit_HEC_constp, " \n len ", len(Rcrit_HEC_constp))
print(" \n Rcrit_IC ", Rcrit_IC_constp, " \n len ", len(Rcrit_IC_constp))

Rcrit_constp_min =  []

for ii in range(0, 37, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("constp", ii+1)

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[0]))

  # print(" \n t_arr ",t_arr)

  Rcrit_constp = h.critical_radius(t_arr, p_arr[0])

  # print(" \n Rcrit_constp ", Rcrit_constp, " \n bool matrix ", ~np.isnan(Rcrit_constp))

  Rcrit_constp_min.append(np.min(Rcrit_constp[~np.isnan(Rcrit_constp)]))

  # print(" \n min Rcrit ", Rcrit_constp[~np.isnan(Rcrit_constp)])

  # Rcrit_HEC = h.critical_radius(t_arr, p_arr)
  # Rcrit_IC = h.critical_radius(t_IC_constp, p_arr)


Rcrit_constp_min = np.array(Rcrit_constp_min)

print(" \nRcrit_min looks like ", Rcrit_constp_min, " \n len ", len(Rcrit_constp_min))

##################################################################
##                      cosnt-Q data                            ##
##################################################################

p_arr_HEC_constQ = np.array([22.0626, 21.5357, 21.696, 20.919])
p_arr_IC_constQ = np.array([22.0626, 21.5357, 21.696, 20.7277])

T_arr_HEC_constQ = np.array([2.058, 2.102, 2.084, 2.121])
t_arr_HEC_constQ = T_arr_HEC_constQ/(1000*SC_beta.Tcp(p_arr_HEC_constQ))

T_arr_IC_constQ = np.array([2.058, 2.102, 2.084, 2.113])
t_arr_IC_constQ = T_arr_IC_constQ/(1000*SC_beta.Tcp(p_arr_IC_constQ))

Rcrit_HEC_constQ = h.critical_radius(t_arr_HEC_constQ, p_arr_HEC_constQ)
Rcrit_IC_constQ = h.critical_radius(t_arr_IC_constQ, p_arr_IC_constQ)


Rcrit_constQ_min = []

for ii in range(0, 4, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("constQ", ii+1)

  print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  print(" \n t_arr ",t_arr)

  Rcrit_constQ = h.critical_radius(t_arr, p_arr)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_constQ_min.append(np.min(Rcrit_constQ[~np.isnan(Rcrit_constQ)]))

Rcrit_constQ_min = np.array(Rcrit_constQ_min)

print(" \n Rcrit_constQ_min ", Rcrit_constQ_min)

print(" \n Rcrit_HEC_constQ ", Rcrit_HEC_constQ)
print(" \n Rcrit_IC_constQ ", Rcrit_IC_constQ)

##################################################################
##                      others-1 data                           ##
##################################################################

p_arr_O1 = np.array([23., 25., 21., 24.865, 27., 27., 24.])

T_arr_HEC_O1 = np.array([2.012, 1.949, 2.123, 1.951, 1.901, 1.903, 1.959])
t_arr_HEC_O1 = T_arr_HEC_O1/(1000*SC_beta.Tcp(p_arr_O1))

T_arr_IC_O1 = T_arr_HEC_O1
t_arr_IC_O1 = T_arr_IC_O1/(1000*SC_beta.Tcp(p_arr_O1))

Rcrit_HEC_O1 = h.critical_radius(t_arr_HEC_O1, p_arr_O1)
Rcrit_IC_O1 = h.critical_radius(t_arr_IC_O1, p_arr_O1)


Rcrit_O1_min = []

for ii in range(0, 7, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("others_1", ii+1)

  print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  print(" \n t_arr ",t_arr)

  Rcrit_O1 = h.critical_radius(t_arr, p_arr)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_O1_min.append(np.min(Rcrit_O1[~np.isnan(Rcrit_O1)]))

Rcrit_O1_min = np.array(Rcrit_O1_min)

print(" \n Rcrit_constQ_min ", Rcrit_O1_min)

print(" \n Rcrit_HEC_constQ ", Rcrit_HEC_O1)
print(" \n Rcrit_IC_constQ ", Rcrit_IC_O1)

##################################################################
##                      others-2 data                           ##
##################################################################

p_O2 = 23.

T_arr_HEC_O2 = np.array([2.0106, 2.0617, 2.0996])
t_arr_HEC_O2 = T_arr_HEC_O2/(1000*SC_beta.Tcp(p_O2))

T_arr_IC_O2 = T_arr_HEC_O2
t_arr_IC_O2 = T_arr_IC_O2/(1000*SC_beta.Tcp(p_O2))

Rcrit_HEC_O2 = h.critical_radius(t_arr_HEC_O2, p_O2)
Rcrit_IC_O2 = h.critical_radius(t_arr_IC_O2, p_O2)

print(" \n Rcrit_HEC_constQ ", Rcrit_HEC_O2)
print(" \n Rcrit_IC_constQ ", Rcrit_IC_O2)

Rcrit_O2_min = []

for ii in range(0, 2, 1):

  print(" \n ii+1 ", ii+1)
  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("others_2", ii+1)

  print(" \n p_arr ",p_arr," \n T_arr ", T_arr  )

  print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  print(" \n t_arr ",t_arr)

  Rcrit_O2 = h.critical_radius(t_arr, p_arr)

  print(" \n Rcrit_O2 ", Rcrit_O2)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_O2_min.append(np.min(Rcrit_O2[~np.isnan(Rcrit_O2)]))

Rcrit_O2_min = np.array(Rcrit_O2_min)

print(" \n Rcrit_O2_min ", Rcrit_O2_min)


  
##################################################################
##                        data plot                             ##
##################################################################

fig, ax = plt.subplots(1, 1)
mSize = 90
unitC = 10**(-3) # nm to micron

s1 = ax.scatter(1/(unitC*Rcrit_constp_min[0:(37-3)]), 1/(unitC*Rcrit_HEC_constp[0:(39-5)]), color = "green", marker = "s", s = mSize+10, label="Const-p, HEC")
s2 = ax.scatter(1/(unitC*Rcrit_constp_min[0:(37-3)]), 1/(unitC*Rcrit_IC_constp[0:(39-5)]), color = "orange", marker = "o", s = mSize, label="Const-p, IC")

s3 = ax.scatter(1/(unitC*Rcrit_constQ_min), 1/(unitC*Rcrit_IC_constQ), color = "red", marker = "h", s = mSize+40, label="Const-Q, IC")
s4 = ax.scatter(1/(unitC*Rcrit_constQ_min), 1/(unitC*Rcrit_HEC_constQ), color = "c", marker = "P", s = mSize, label="Const-Q, HEC")

s5 = ax.scatter(1/(unitC*Rcrit_O1_min), 1/(unitC*Rcrit_IC_O1), color = "turquoise", marker = "v", s = mSize+20, label="other-1, IC")
s6 = ax.scatter(1/(unitC*Rcrit_O1_min), 1/(unitC*Rcrit_HEC_O1), color = "violet", marker = "*", s = mSize, label="other-1, HEC")

print(" \n Rcrit_IC_O1[0:2] ", Rcrit_IC_O1[0:2])
s7 = ax.scatter(1/(unitC*Rcrit_O2_min), 1/(unitC*Rcrit_IC_O1[0:2]), color = "plum", marker = "d", s = mSize+20, label="other-2, IC")
s8 = ax.scatter(1/(unitC*Rcrit_O2_min), 1/(unitC*Rcrit_HEC_O1[0:2]), color = "indigo", marker = "+", s = mSize, label="other-2, HEC")


ax.set_xlabel(r"$R^{-1(A)}_{crit}/{\mu}m^{-1}$", fontsize=20.0)
ax.set_ylabel(r"$R^{-1(B)}_{crit}/{\mu}m^{-1}$", fontsize=20.0)
# ax.set_xlim([0.0, 0.6])
# ax.set_ylim([0.0, 1.5])
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)
ax.grid(True)
leg = ax.legend([s1, s2, s3, s4, s5, s6, s7, s8],["Const-p, HEC", "Const-p, IC", "Const-Q, IC", "Const-Q, HEC", "Others-1, IC","Others-1, HEC", "Others-2, IC","Others-2, HEC"], fontsize=15.0, loc='lower right')

plt.show()
