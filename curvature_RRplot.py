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
import he3_tools_Vn01 as hn

import Module_RR_pT as rr_pT
import Module_SC_Beta_V05 as SC_beta

import Module_New_Parpia_pT_parameter as Parpia

h.set_default("DEFAULT_SC_ADJUST", True)
h.set_default("DEFAULT_T_SCALE", "PLTS")

SC_beta.turn_on_PLTS()

################################################################
##                       cosnt-p data                         ##
################################################################


p_arr_constp = Parpia.pressure
T_IC_constp = Parpia.T_IC; t_IC_constp = T_IC_constp/(1000*SC_beta.Tcp(p_arr_constp))
T_HEC_constp =Parpia.T_HEC; t_HEC_constp = T_HEC_constp/(1000*SC_beta.Tcp(p_arr_constp))

# print(" \n t_IC_constp looks like ", t_IC_constp)

# Rcrit in unit of micron
Rcrit_HEC_constp = (10**(-3))*h.critical_radius(t_HEC_constp, p_arr_constp)
Rcrit_IC_constp = (10**(-3))*h.critical_radius(t_IC_constp, p_arr_constp)

# in unit of micron
xiGL_IC_constp = (10**6)*(SC_beta.xiGL_JWS(p_arr_constp, T_IC_constp*(10**(-3))))
xiGL_HEC_constp =(10**6)*(SC_beta.xiGL_JWS(p_arr_constp, T_HEC_constp*(10**(-3))))

print(" \n Rcrit_HEC_constp ", Rcrit_HEC_constp, " \n len ", len(Rcrit_HEC_constp))
print(" \n Rcrit_IC_constp ", Rcrit_IC_constp, " \n len ", len(Rcrit_IC_constp))

Rcrit_constp_min = []
xiGL_constp_min = []

for ii in range(0, 39, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("constp", ii+1)

  # print(" \n p_arr looks like ",p_arr, " \n T_arr looks like ", T_arr)

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[0]))

  # print(" \n t_arr ",t_arr)

  Rcrit_constp = (10**(-3))*h.critical_radius(t_arr, p_arr[0])

  # print(" \n Rcrit_constp ", Rcrit_constp, " \n bool matrix ", ~np.isnan(Rcrit_constp))

  Rcrit_constp_min.append(np.min(Rcrit_constp[~np.isnan(Rcrit_constp)]))
  
  boolen_constp = Rcrit_constp[~np.isnan(Rcrit_constp)] == np.min(Rcrit_constp[~np.isnan(Rcrit_constp)])
  print(" \n boolen_constp looks like ", boolen_constp)

  # print(" \n xiGL arr looks like ", (10**6)*SC_beta.xiGL_JWS(p_arr, T_arr*(10**(-3))))
  
  xiGL_constp_min.append((((10**6)*SC_beta.xiGL_JWS(p_arr, T_arr*(10**(-3))))[~np.isnan(Rcrit_constp)])[boolen_constp][0])

  # print(" \n min Rcrit ", Rcrit_constp[~np.isnan(Rcrit_constp)])

  # Rcrit_HEC = h.critical_radius(t_arr, p_arr)
  # Rcrit_IC = h.critical_radius(t_IC_constp, p_arr)


Rcrit_constp_min = np.array(Rcrit_constp_min)
xiGL_constp_min = np.array(xiGL_constp_min)

print(" \nRcrit_min looks like ", Rcrit_constp_min, " \n len ", len(Rcrit_constp_min))
print(" \n xiGL_constp__min ", xiGL_constp_min)


##################################################################
##                      cosnt-Q data                            ##
##################################################################

# all must be PLTS scale
p_arr_HEC_constQ = np.array([22.0626, 21.5357, 21.696, 20.919])
p_arr_IC_constQ = np.array([22.0626, 21.5357, 21.696, 20.7277])

T_arr_HEC_constQ = np.array([2.058, 2.102, 2.084, 2.121])
t_arr_HEC_constQ = T_arr_HEC_constQ/(1000*SC_beta.Tcp(p_arr_HEC_constQ))

T_arr_IC_constQ = np.array([2.058, 2.102, 2.084, 2.113])
t_arr_IC_constQ = T_arr_IC_constQ/(1000*SC_beta.Tcp(p_arr_IC_constQ))

# Rcrit in unit of micron
Rcrit_HEC_constQ = (10**(-3))*h.critical_radius(t_arr_HEC_constQ, p_arr_HEC_constQ)
Rcrit_IC_constQ = (10**(-3))*h.critical_radius(t_arr_IC_constQ, p_arr_IC_constQ)

# xiGL
xiGL_HEC_constQ = (10**6)*SC_beta.xiGL_JWS(p_arr_HEC_constQ, T_arr_HEC_constQ*(10**(-3)))
xiGL_IC_constQ = (10**6)*SC_beta.xiGL_JWS(p_arr_IC_constQ, T_arr_IC_constQ*(10**(-3)))


Rcrit_constQ_min = []
xiGL_constQ_min = []

for ii in range(0, 4, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("constQ", ii+1)

  # print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  # print(" \n t_arr ",t_arr)

  Rcrit_constQ = (10**(-3))*h.critical_radius(t_arr, p_arr)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_constQ_min.append(np.min(Rcrit_constQ[~np.isnan(Rcrit_constQ)]))

  boolen_constQ = Rcrit_constQ[~np.isnan(Rcrit_constQ)] == np.min(Rcrit_constQ[~np.isnan(Rcrit_constQ)])

  xiGL_constQ_min.append((((10**6)*SC_beta.xiGL_JWS(p_arr, T_arr*(10**(-3))))[~np.isnan(Rcrit_constQ)])[boolen_constQ][0])

Rcrit_constQ_min = np.array(Rcrit_constQ_min)
xiGL_constQ_min = np.array(xiGL_constQ_min)

print(" \n Rcrit_constQ_min ", Rcrit_constQ_min)
print(" \n xiGL_constQ_min ", xiGL_constQ_min)

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

Rcrit_HEC_O1 = (10**(-3))*h.critical_radius(t_arr_HEC_O1, p_arr_O1)
Rcrit_IC_O1 = (10**(-3))*h.critical_radius(t_arr_IC_O1, p_arr_O1)

# xiGL
xiGL_HEC_O1 = (10**6)*SC_beta.xiGL_JWS(p_arr_O1, T_arr_HEC_O1*(10**(-3)))
xiGL_IC_O1 = (10**6)*SC_beta.xiGL_JWS(p_arr_O1, T_arr_IC_O1*(10**(-3)))



Rcrit_O1_min = []
xiGL_O1_min = []

for ii in range(0, 7, 1):

  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("others_1", ii+1)

  # print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  # print(" \n t_arr ",t_arr)

  Rcrit_O1 = (10**(-3))*h.critical_radius(t_arr, p_arr)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_O1_min.append(np.min(Rcrit_O1[~np.isnan(Rcrit_O1)]))

  boolen_O1 = Rcrit_O1[~np.isnan(Rcrit_O1)] == np.min(Rcrit_O1[~np.isnan(Rcrit_O1)])

  xiGL_O1_min.append((((10**6)*SC_beta.xiGL_JWS(p_arr, T_arr*(10**(-3))))[~np.isnan(Rcrit_O1)])[boolen_O1][0])

Rcrit_O1_min = np.array(Rcrit_O1_min)
xiGL_O1_min = np.array(xiGL_O1_min)

print(" \n Rcrit_other-1_min ", Rcrit_O1_min)
print(" \n xiGL_O1_min ", xiGL_O1_min)

print(" \n Rcrit_HEC_other-1 ", Rcrit_HEC_O1)
print(" \n Rcrit_IC_other-1 ", Rcrit_IC_O1)

##################################################################
##                      others-2 data                           ##
##################################################################

p_O2 = 23.

T_arr_HEC_O2 = np.array([2.0106, 2.0617, 2.0996])
t_arr_HEC_O2 = T_arr_HEC_O2/(1000*SC_beta.Tcp(p_O2))

T_arr_IC_O2 = T_arr_HEC_O2
t_arr_IC_O2 = T_arr_IC_O2/(1000*SC_beta.Tcp(p_O2))

Rcrit_HEC_O2 = (10**(-3))*h.critical_radius(t_arr_HEC_O2, p_O2)
Rcrit_IC_O2 = (10**(-3))*h.critical_radius(t_arr_IC_O2, p_O2)

# xiGL
xiGL_HEC_O2 = (10**6)*SC_beta.xiGL_JWS(p_O2, T_arr_HEC_O2*(10**(-3)))
xiGL_IC_O2 = (10**6)*SC_beta.xiGL_JWS(p_O2, T_arr_IC_O2*(10**(-3)))

print(" \n Rcrit_HEC_other-2 ", Rcrit_HEC_O2)
print(" \n Rcrit_IC_other-2 ", Rcrit_IC_O2)

print("\n xiGL_HEC_O2 ", xiGL_HEC_O2)
print("\n xiGL_IC_O2 ", xiGL_IC_O2)

Rcrit_O2_min = []
xiGL_O2_min = []

for ii in range(0, 3, 1):

  print(" \n ii+1 looks like ", ii+1)
  # load the rr_pT data  
  p_arr, T_arr = rr_pT.get_data("others_2", ii+1)

  # print(" \n p_arr ",p_arr," \n T_arr ", T_arr  )

  # print(" \n p_arr[-1] ", p_arr[-1])

  t_arr = T_arr/(1000*SC_beta.Tcp(p_arr[-1]))

  # print(" \n t_arr ",t_arr)

  Rcrit_O2 = (10**(-3))*h.critical_radius(t_arr, p_arr)

  # print(" \n Rcrit_O2 ", Rcrit_O2)

  # print(" \n Rcrit_constQ ", Rcrit_constQ, " \n bool matrix ", ~np.isnan(Rcrit_constQ))

  Rcrit_O2_min.append(np.min(Rcrit_O2[~np.isnan(Rcrit_O2)]))

  boolen_O2 = Rcrit_O2[~np.isnan(Rcrit_O2)] == np.min(Rcrit_O2[~np.isnan(Rcrit_O2)])

  xiGL_O2_min.append((((10**6)*SC_beta.xiGL_JWS(p_arr, T_arr*(10**(-3))))[~np.isnan(Rcrit_O2)])[boolen_O2][0])

Rcrit_O2_min = np.array(Rcrit_O2_min)
xiGL_O2_min = np.array(xiGL_O2_min)

print(" \n Rcrit_O2_min ", Rcrit_O2_min)
print(" \n xiGL_O2_min ", xiGL_O2_min)


  
##################################################################
##                        data plot                             ##
##################################################################

##################################################################
##        >>>>>>>>>>>>>>  curvature plot  <<<<<<<<<<<<<         ##
##################################################################

fig, ax = plt.subplots(1, 1)
mSize = 90


s1 = ax.scatter(1/(Rcrit_constp_min), 1/(Rcrit_HEC_constp), color = "green", marker = "s", s = mSize+10, label="Const-p, HEC")
s2 = ax.scatter(1/(Rcrit_constp_min), 1/(Rcrit_IC_constp), color = "orange", marker = "o", s = mSize, label="Const-p, IC")

# axt = ax.twinx()
# g1 = axt.scatter(1/(Rcrit_constp_min[0:(37-3)]), xiGL_HEC_constp[0:(37-3)], color = 'red', marker = 'o', s=10)
# g2 = axt.scatter(1/(Rcrit_constp_min[0:(37-3)]), xiGL_IC_constp[0:(37-3)], color = 'blue', marker = 's', s=10)

s3 = ax.scatter(1/(Rcrit_constQ_min), 1/(Rcrit_IC_constQ), color = "red", marker = "h", s = mSize+80, label="Const-Q, IC")
s4 = ax.scatter(1/(Rcrit_constQ_min), 1/(Rcrit_HEC_constQ), color = "c", marker = "P", s = mSize, label="Const-Q, HEC")

s5 = ax.scatter(1/(Rcrit_O1_min), 1/(Rcrit_IC_O1), color = "turquoise", marker = "v", s = mSize+20, label="other-1, IC")
s6 = ax.scatter(1/(Rcrit_O1_min), 1/(Rcrit_HEC_O1), color = "violet", marker = "*", s = mSize, label="other-1, HEC")

# print(" \n Rcrit_IC_O1[0:2] ", Rcrit_IC_O1[0:2])
s7 = ax.scatter(1/(Rcrit_O2_min), 1/(Rcrit_IC_O2), color = "plum", marker = "d", s = mSize+20, label="other-2, IC")
s8 = ax.scatter(1/(Rcrit_O2_min), 1/(Rcrit_HEC_O2), color = "indigo", marker = "4", s = mSize+80, label="other-2, HEC")


ax.set_xlabel(r"${\kappa}^{(A)}_{crit}/{\mu}m^{-1}$", fontsize=20.0)
ax.set_ylabel(r"${\kappa}^{(B)}_{crit}/{\mu}m^{-1}$", fontsize=20.0)
ax.set_xlim([0.0, 0.7])
ax.set_ylim([0.0, 1.4])

# axt.set_ylabel(r"$\xi_{GL}(p, T)/{\mu}m$")

plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)
ax.grid(True)
# leg = ax.legend([s1, s2, s3, s4, s5, s6, s7, s8],["Const-p, HEC", "Const-p, IC", "Const-Q, IC", "Const-Q, HEC", "Others-1, IC","Others-1, HEC", "Others-2, IC","Others-2, HEC"], fontsize=15.0, loc='lower right')
leg = ax.legend([s1, s2, s3, s4, s5, s6, s7, s8],["Const-p, HEC", "Const-p, IC", "Const-Q, IC", "Const-Q, HEC", "Others-1, IC","Others-1, HEC", "Other-2, IC", "Other-2, HEC"], fontsize=15.0, loc='upper left')

# print(" \n kappa * xiGL: IC \n ", xiGL_IC_constp[0:(37-3)]*(1/(unitC*Rcrit_constp_min[0:(37-3)])))

plt.show()


####################################################################
##                   \kappa \xiGl vs kappa \xiGL                  ##
####################################################################

fig4, ax4 = plt.subplots(1,1)

KA_cp = 1/(Rcrit_constp_min)
KB_IC_cp = 1/(Rcrit_IC_constp)
KB_HEC_cp = 1/(Rcrit_HEC_constp)

KA_cQ = 1/(Rcrit_constQ_min)
KB_IC_cQ = 1/(Rcrit_IC_constQ)
KB_HEC_cQ = 1/(Rcrit_HEC_constQ)

KA_O1 = 1/(Rcrit_O1_min)
KB_IC_O1 = 1/(Rcrit_IC_O1)
KB_HEC_O1 = 1/(Rcrit_HEC_O1)

KA_O2 = 1/(Rcrit_O2_min)
KB_IC_O2 = 1/(Rcrit_IC_O2)
KB_HEC_O2 = 1/(Rcrit_HEC_O2)

KAxiGL_cp = KA_cp*(xiGL_constp_min)
print("\n KAxiGL: ",KAxiGL_cp, " \n KAxiGL.shape ",KAxiGL_cp.shape)

KBICxiGL_cp = KB_IC_cp*(xiGL_IC_constp)
print("\n KBICxiGL: ",KBICxiGL_cp, " \n KBICxiGL.shape ",KBICxiGL_cp.shape)

KBHECxiGL_cp = KB_HEC_cp*(xiGL_HEC_constp)
print("\n KBHECxiGL: ",KBICxiGL_cp, " \n KBHECxiGL.shape ",KBHECxiGL_cp.shape)


KAxiGL_cQ = KA_cQ*(xiGL_constQ_min)
print("\n KAxiGL: ",KAxiGL_cQ, " \n KAxiGL.shape ",KAxiGL_cQ.shape)

KBICxiGL_cQ = KB_IC_cQ*(xiGL_IC_constQ)
print("\n KBICxiGL: ",KBICxiGL_cQ, " \n KBICxiGL.shape ",KBICxiGL_cQ.shape)

KBHECxiGL_cQ = KB_HEC_cQ*(xiGL_HEC_constQ)
print("\n KBHECxiGL: ",KBICxiGL_cQ, " \n KBHECxiGL.shape ",KBHECxiGL_cQ.shape)


KAxiGL_O1 = KA_O1*(xiGL_O1_min)
print("\n KAxiGL: ",KAxiGL_O1, " \n KAxiGL.shape ",KAxiGL_O1.shape)

KBICxiGL_O1 = KB_IC_O1*(xiGL_IC_O1)
print("\n KBICxiGL: ",KBICxiGL_O1, " \n KBICxiGL.shape ",KBICxiGL_O1.shape)

KBHECxiGL_O1 = KB_HEC_O1*(xiGL_HEC_O1)
print("\n KBHECxiGL: ",KBICxiGL_O1, " \n KBHECxiGL.shape ",KBHECxiGL_O1.shape)


KAxiGL_O2 = KA_O2*(xiGL_O2_min)
print("\n KAxiGL: ",KAxiGL_O2, " \n KAxiGL.shape ",KAxiGL_O2.shape)

KBICxiGL_O2 = KB_IC_O2*(xiGL_IC_O2)
print("\n KBICxiGL: ",KBICxiGL_O2, " \n KBICxiGL.shape ",KBICxiGL_O2.shape)

KBHECxiGL_O2 = KB_HEC_O2*(xiGL_HEC_O2)
print("\n KBHECxiGL: ",KBICxiGL_O2, " \n KBHECxiGL.shape ",KBHECxiGL_O2.shape)


# print("\n xiGL_constp_min.shape: ", xiGL_constp_min[0:(37-3)].shape, " KA.shape: ", KA.shape, " KB_IC.shape: ", KB_IC.shape)

q1 = ax4.scatter(KAxiGL_cp, KBICxiGL_cp, color = "green", marker = "s", s = mSize+10, label="Const-p, IC")
q2 = ax4.scatter(KAxiGL_cp, KBHECxiGL_cp, color = "darkorange", marker = "o", s = mSize, label="Const-p, HEC")

q3 = ax4.scatter(KAxiGL_cQ, KBICxiGL_cQ, color = "blue", marker = "p", s = mSize+60, label="Const-Q, IC")
q4 = ax4.scatter(KAxiGL_cQ, KBHECxiGL_cQ, color = "deeppink", marker = "h", s = mSize, label="Const-Q, HEC")

q5 = ax4.scatter(KAxiGL_O1, KBICxiGL_O1, color = "chocolate", marker = "*", s = mSize+60, label="Const-Q, IC")
q6 = ax4.scatter(KAxiGL_O1, KBHECxiGL_O1, color = "lime", marker = "x", s = mSize, label="Const-Q, HEC")

q7 = ax4.scatter(KAxiGL_O2, KBICxiGL_O2, color = "purple", marker = "<", s = mSize+60, label="Const-Q, IC")
q8 = ax4.scatter(KAxiGL_O2, KBHECxiGL_O2, color = "fuchsia", marker = "4", s = mSize+60, label="Const-Q, HEC")

ax4.set_xlabel(r"$\kappa^{(A)}_{crit}{\xi_{GL}(p,T)}$", fontsize=20.0)
ax4.set_ylabel(r"$\kappa^{(B)}_{crit}{\xi_{GL}(p,T)}$", fontsize=20.0)
ax4.set_xlim([0.0, 0.04])
ax4.set_ylim([0.0, 0.03])
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)

ax4.grid(True)

leg4 = ax4.legend([q1, q2, q3, q4, q5, q6, q7, q8],["Const-p, HEC", "Const-p, IC", "Const-Q, IC", "Const-Q, HEC", "Others-1, IC","Others-1, HEC", "Other-2, IC", "Other-2, HEC"], fontsize=15.0, loc='upper right')


plt.show()

####################################################################
##             >>>>>>>>>>>>>>  radius plot  <<<<<<<<<<<<          ##
####################################################################

fig1, ax1 = plt.subplots(1,1)

mSize = 90


ss1 = ax1.scatter(Rcrit_constp_min[0:(37-3)], Rcrit_HEC_constp[0:(39-5)], color = "green", marker = "s", s = mSize+10, label="Const-p, HEC")
ss2 = ax1.scatter(Rcrit_constp_min[0:(37-3)], Rcrit_IC_constp[0:(39-5)], color = "orange", marker = "o", s = mSize, label="Const-p, IC")

print(" \n Rcrit_IC_constQ = ", Rcrit_IC_constQ, " \n Rcrit_HEC_constQ ", Rcrit_HEC_constQ)
ss3 = ax1.scatter(Rcrit_constQ_min, Rcrit_IC_constQ, color = "red", marker = "h", s = mSize+50, label="Const-Q, IC")
ss4 = ax1.scatter(Rcrit_constQ_min, Rcrit_HEC_constQ, color = "c", marker = "P", s = mSize, label="Const-Q, HEC")

ss5 = ax1.scatter(Rcrit_O1_min, Rcrit_IC_O1, color = "turquoise", marker = "v", s = mSize+20, label="other-1, IC")
ss6 = ax1.scatter(Rcrit_O1_min, Rcrit_HEC_O1, color = "violet", marker = "*", s = mSize, label="other-1, HEC")

# # print(" \n Rcrit_IC_O1[0:2] ", Rcrit_IC_O1[0:2])
ss7 = ax1.scatter(Rcrit_O2_min, Rcrit_IC_O2, color = "plum", marker = "d", s = mSize+20, label="other-2, IC")
ss8 = ax1.scatter(Rcrit_O2_min, Rcrit_HEC_O2, color = "indigo", marker = "+", s = mSize, label="other-2, HEC")


ax1.set_xlabel(r"$R^{(A)}_{crit}/{\mu}m$", fontsize=20.0)
ax1.set_ylabel(r"$R^{(B)}_{crit}/{\mu}m$", fontsize=20.0)
ax1.set_xlim([0.0, 25])
ax1.set_ylim([0.0, 14])
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)
ax1.grid(True)

leg1 = ax1.legend([ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8],["Const-p, HEC", "Const-p, IC", "Const-Q, IC", "Const-Q, HEC", "Others-1, IC","Others-1, HEC", "Others-2, IC","Others-2, HEC"], fontsize=15.0, loc='upper left')


plt.show()
