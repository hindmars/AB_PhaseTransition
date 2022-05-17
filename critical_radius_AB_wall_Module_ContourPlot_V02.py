############################################################################
######               Important notations and descriptions            #######
############################################################################

#           >>>>>>>>>>>>>    Read Me, please <<<<<<<<<<<<<<<<

# This script is used for calculating and plotting the thin-wall evaluation of critical radius half-cylinder bubble,
# which is supposed appears in Parpia group's experiments when B-phase in HEC chamber propagate to IC chamber.

# This script use *SC_beta, *AB-wall and *he3_tools modules. To solve the discrapcy of T_AB line between RWS SC and Greywall data,
# a pressure-dependent temperature shift is added, then the superheating results are unrelaible. 

# This is version 0.1 

# author: Quang (timohyva@github)

############################################################################
######                       Modules & Constants                     #######
############################################################################

import Module_SC_Beta_V05 as SCC
import Module_AB_wall_V00 as WB
import he3_tools_Vn01 as h

# Tc, TAB data digitized from Parpia's manuscript
# import Module_p_Tc_TAB_Parpia as pTcTAB

# import Module_plot_TAB_line as TAB_RWS
import Module_Lotynk_pT_parameter as Lotynk
import Module_New_Parpia_pT_parameter as Parpia
import Module_Parpia_pT_ConsQ as ConsQ

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from math import pi
import math

# main code starts from here

# Length unit, Time unit, Energy unit, mass unit, pressure unit 
zeta3 = 1.2020569;
m = 1;s = 1; J = 1; Kelvin = 1; kg =1; bar = 1

kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1

hbar = 1.054571817*(10**(-34)) # planck constant, J.s

u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom
              
##############################################################################
#####                    calculate the Rc array                          #####
##############################################################################

Pressure = np.arange(19.0,30.1,0.1) # bar
# Pressure = np.arange(20.0,30.1,0.1) # bar
# Pressure = np.arange(20.0, 21.2, 0.2)

Temperature = (np.arange(1.50, 2.5, 0.01))*(10**(-3)) # Kelvin

# TAB_RWSco = h.T_mK(h.t_AB(Pressure), Pressure)

# TAB_Greywall = h.TAB_poly_Greywall(Pressure)

#####################################################################
###           mask the TAB_greywall poly when p < p_pcp          ####
#####################################################################

TAB_Greywall_arr = np.array([])

for p in Pressure:
  dp = p - h.p_pcp_bar

  if dp < 0.:
    TAB_Greywall_arr = np.append(TAB_Greywall_arr,np.nan)
  elif dp >= 0.:
     TAB_Greywall_arr = np.append(TAB_Greywall_arr, h.TAB_poly_Greywall(dp))
     
#####################################################################

#####################################################################
###           mask the TAB_greywall poly when p < p_pcp          ####
#####################################################################

TAB_PLTS_arr = np.array([])

for p in Pressure:
  # (h.TAB_poly_PLTS(Pressure)

  if p < h.p_pcp_bar:
    TAB_PLTS_arr = np.append(TAB_PLTS_arr,np.nan)
  elif p >= h.p_pcp_bar:
    TAB_PLTS_arr = np.append(TAB_PLTS_arr, h.TAB_poly_PLTS(p))
     
######################################################################

T_array, P_array = np.meshgrid(Temperature, Pressure)

print("\n P_array looks like\n", P_array,"\n \n T_array looks like \n", T_array)

# print("\n TAB_RWS looks like\n", TAB_RWS,"\n TAB_Greywall looks like\n", TAB_Greywall)

print("\n TAB_Greywall looks like\n", TAB_Greywall_arr)



##############################################################################
#####                   calcuate sigmaAB for pTcTAB.p                    #####
##############################################################################

tensionAB_p = np.array([])
for p in Pressure:
  if p < 21.5:
   tensionAB_p = np.append(tensionAB_p, 0.734)
     
  elif p >= 21.5: 
   tensionAB_p = np.append(tensionAB_p, WB.get_wall_n_surfaceEnergy(p))

print("\n tensionAB_p looks like ", tensionAB_p)


##############################################################################


Rc = np.zeros((len(Pressure),len(Temperature)))
Rc_sh = np.zeros(Rc.shape)


# #############################################################################

for iP in range(0, len(Pressure), 1):
    print('\n\n Now P is:', Pressure[iP], '\n\n')

    p = Pressure[iP]


    sigma_p = tensionAB_p[iP]

    for it in range(0, len(Temperature), 1):

       T = Temperature[it]

       t = T/SCC.Tcp(p)

       print('\n now temperature is:, ', T)

       if (p < 21.22) and (t <= 1.):

        fB = SCC.alpha_td(p,T)/SCC.betaB_td(p,T)
     
        fBA = (-SCC.alpha_td(p,T)/SCC.betaB_td(p,T)) + (SCC.alpha_td(p,T)/SCC.betaA_td(p,T))
        
        xiGL = SCC.xiGL_OC(p, T)

        Rc[iP, it] =  abs(-1.*((sigma_p*fB)/fBA)*xiGL)
        Rc_sh[iP, it] = np.nan

        print(" \n pressure lower than P_PCP, sigma/fBA looks like, ", Rc[iP, it]*(10**(6)))

       else: 
    
         if (t >= 1.):

           print(" bro, we just got temperature at Tc, save a np.nan. ")

           # mask data from t > 1 region

           Rc[iP, it] = np.nan

           Rc_sh[iP, it] = np.nan

      # only for A to B  
      # elif (t <1.) and (t >= h.t_AB(p)):

      # A to B or B to A
         elif (t == SCC.tAB_RWSco(p)):

            print(" bro, you got on T_AB, radius divergence. ")

            Rc[iP, it] = np.nan
            Rc_sh[iP, it] = np.nan
          
                 
         else:
        
            fB = SCC.alpha_td(p,T)/SCC.betaB_td(p,T)
     
            fBA = (-SCC.alpha_td(p,T)/SCC.betaB_td(p,T)) + (SCC.alpha_td(p,T)/SCC.betaA_td(p,T))
        
            xiGL = SCC.xiGL_OC(p, T)

            print(" \n fB is ", fB, " fBA is ", fBA, " xiGL is ", xiGL, " sigmaAB is", sigma_p)

            if t < SCC.tAB_RWSco(p):

              Rc[iP, it] = abs(-1.*((sigma_p*fB)/fBA)*xiGL)
              Rc_sh[iP, it] = np.nan
              print("\n t < t_AB, sigmaAB/fBA looks like : ", Rc[iP, it])
          
            elif t > SCC.tAB_RWSco(p):

              Rc_sh[iP, it] = abs(-1.*((sigma_p*fB)/fBA)*xiGL)
              Rc[iP, it] = np.nan        
          

# -1sigma/fBA
# print("\n tensionAB looks like",tensionAB)
# Rc = -1.*((tensionAB*fB)/fBA)*xiGL; print("\n sigmaAB/fBA looks like : ", Rc)

print(" \n Rc array looks like ",Rc, "\n Rc_sh array looks like ", Rc_sh)




#############################################################################
#####           scatter plot the ratio verse the pressure            ########
#############################################################################.

fig, ax = plt.subplots(1,1)
Levels1 = [0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2., 2.5, 3., 4., 5., 6.]
Levels2 = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 2., 2.5, 3, 4., 5., 6., 7., 8., 9., 10., 11]

Levels1sh = [0.1, 0.3, 0.5, 0.9, 1., 1.2, 1.4, 1.7, 2.0]

# plot the TAB line with RWS2019 SC
# TABrws = ax.contour(TAB_RWS.X*(10**3), TAB_RWS.Y, TAB_RWS.EnergyDensity_Difference_fABGL, levels=[0.0], colors='red')

cf1 = ax.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, (1/4)*Rc*(10**(6)), cmap=cm.PuBu_r, levels=Levels2)
ax.set_ylabel(r'$p/bar$', fontsize=15.0,);
ax.set_xlabel(r'$T$/mK', fontsize=15.0,);
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)

fig.colorbar(cf1, ax=ax, location = 'left')

# cf2 = ax.contourf(h.T_GtoPlts6_poly(T_array*(10**(3)) + T_shift_array), P_array, Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels2)
cf2 = ax.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, (1/4)*Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels1sh)

fig.colorbar(cf2, ax=ax, location = 'right')

c3 = ax.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, (1/4)*Rc*(10**(6)), levels=Levels1, colors='orange')
plt.clabel(c3, inline=True, fontsize=10.5, colors='r')

# c4 = ax.contour(h.T_GtoPlts6_poly(T_array*(10**(3)) +  T_shift_array), P_array, Rc_sh*(10**(6)), levels=Levels1, colors='blue')
c4 = ax.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, (1/4)*Rc_sh*(10**(6)), levels=Levels1sh, colors='blue')
plt.clabel(c4, inline=True, fontsize=10.5, colors='b')

# plot Greywall T_AB line
# TABGreywall, = ax.plot(h.T_AB_Greywall_poly(Pressure), Pressure, color = "orange")

# plot PLTS T_AB line from Greywall data
# TABGreywalltoPlts, = ax.plot(h.T_GtoPlts6_poly(h.T_AB_Greywall_poly(Pressure)), Pressure, color = "red")

# plot PLTS T_AB lien
# TABPlts, = ax.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "pink")
TABPlts, = ax.plot(TAB_PLTS_arr, Pressure, color = "orange")

# plot PLTS Tc
# TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
TCPlts, = ax.plot(h.T_GtoPlts6_low_poly(h.Tc_mK(Pressure)), Pressure, color = "purple")

mSize = 140

# scatter plot of parpia's constant pressure data
sIC1 = ax.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label="p-T IC")
sHEC1 = ax.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="p-T HEC")

# scatter plot of parpia's constant Q data

sHEC_CQ1 = ax.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="p-T HEC-CQ")
sIC_CQ1 = ax.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="p-T IC-CQ")

# ax.legend(labels=(r"$T_{AB}$", "p-T IC", "p-T HEC"))

# leg1 = ax.legend(loc="upper right")
# leg2 = ax.legend([TCPlts,TABGreywalltoPlts,sIC,sHEC],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC"],  fontsize=20.5, loc='lower right')
leg2 = ax.legend([TCPlts,TABPlts,sIC1,sHEC1,sHEC_CQ1,sIC_CQ1],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC", r"$p-T-IC, ConsQ$",r"$p-T-HEC, ConsQ$"],  fontsize=15.0, loc='lower right')
# ax.add_artist(leg1)

# ax.set_title(r"contour plot of $(\sigma_{AB}(p)/f_{AB}(p,T))/{\mu}m$")
ax.set_title(r"contour plot of $|(\sigma_{AB}(p)/4f_{AB}(p,T))|/{\mu}m$", fontsize=14.0,)
ax.set_xlim([1.8, 2.4])
ax.set_ylim([19., 30.])
ax.grid(True)

fig.savefig("Rc_contour_RWS_GreywelltoPltsTC_AToBBToA3.pdf")

plt.show()


