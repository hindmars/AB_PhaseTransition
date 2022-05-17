############################################################################
######               Important notations and descriptions            #######
############################################################################

#   >>>>>>>>>>>>>    Read Me, please <<<<<<<<<<<<<<<<

'''This script is used for calculating and plotting the thin-wall evaluation of critical radius half-cylinder bubble,which is supposed appears in Parpia group's experiments when B-phase in HEC chamber propagate to IC chamber.

This script use *SC_beta, *AB-wall and *he3_tools modules. To solve the discrapcy of T_AB line between RWS SC and Greywall data,
a pressure-dependent temperature shift is added, then the superheating results are unrelaible. 

This is version 0.4, after error of alpha_td^2 in fA and fB was found.

author: Quang (timohyva@github)'''

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
zeta3 = 1.2020569031595942854
m = 1;s = 1; J = 1; Kelvin = 1; kg =1; bar = 1

kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1

hbar = 1.054571817*(10**(-34)) # planck constant, J.s

u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom
              
##############################################################################
#####                    calculate the Rc array                          #####
##############################################################################

p_step = 0.5
Pressure = np.arange(20.0,30.0+p_step,0.5) # bar
# Pressure = np.arange(20.0,30.1,0.1) # bar
# Pressure = np.arange(20.0, 21.2, 0.2)

Temperature = (np.arange(1.80, 2.5, 0.01))*(10**(-3)) # Kelvin

# TAB_RWSco = h.T_mK(h.t_AB(Pressure), Pressure)

# TAB_Greywall = h.TAB_poly_Greywall(Pressure)

# IC chamble texture size R
R_ratio = 5. # in unit of \xiGL(p, T) or \xiGL(p)
# R_ratio = 0.0001 # in unit of \xiGL(p, T) or \xiGL(p)

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

# tensionAB_p = np.ones(Pressure.shape)*0.734

##############################################################################


Rc = np.zeros((len(Pressure),len(Temperature)))
Rc_sh = np.zeros(Rc.shape)

RcT = np.zeros((len(Pressure),len(Temperature)))
RcT_sh = np.zeros(Rc.shape)

Nomi_Dino_arr = np.zeros((len(Pressure),len(Temperature)))
Dino_arr = np.zeros((len(Pressure),len(Temperature)))

R_term_arr = np.zeros((len(Pressure),len(Temperature)))
sigma_term_arr = np.zeros((len(Pressure),len(Temperature)))

fBA6_arr = np.zeros((len(Pressure),len(Temperature)))
absfBA_xiGL_arr = np.zeros((len(Pressure),len(Temperature)))


# # #############################################################################

# for iP in range(0, len(Pressure), 1):
#     print('\n\n Now P is:', Pressure[iP], '\n\n')

#     p = Pressure[iP]


#     sigma_p = tensionAB_p[iP]

#     for it in range(0, len(Temperature), 1):

#        T = Temperature[it]

#        t = T/SCC.Tcp(p)

#        print('\n now temperature is:, ', T)

#        if (p < 21.22) and (t <= 1.):
     
#         fBA = SCC.fB(p, T) - SCC.fA(p, T)
        
#         xiGL = SCC.xiGL_OC(p, T)

#         Rc[iP, it] =  -1.*((sigma_p*abs(SCC.fB(p, T))*xiGL)/fBA)
#         Rc_sh[iP, it] = np.nan

#         print(" \n pressure lower than P_PCP, sigma/fBA looks like, ", Rc[iP, it]*(10**(6)))

#        else: 
    
#          if (t >= 1.):

#            print(" bro, we just got temperature at Tc, save a np.nan. ")

#            # mask data from t > 1 region

#            Rc[iP, it] = np.nan

#            Rc_sh[iP, it] = np.nan

#       # only for A to B  
#       # elif (t <1.) and (t >= h.t_AB(p)):

#       # A to B or B to A
#          elif (t == SCC.tAB_RWSco(p)):

#             print(" bro, you got on T_AB, radius divergence. ")

#             Rc[iP, it] = np.nan
#             Rc_sh[iP, it] = np.nan
          
                 
#          else:
        
#             fBA = SCC.fB(p, T) - SCC.fA(p, T)
        
#             xiGL = SCC.xiGL_OC(p, T)


#             print(" \n fB is ", SCC.fB(p, T), " fBA is ", fBA, " xiGL is ", xiGL, " sigmaAB is", sigma_p)

#             if t < SCC.tAB_RWSco(p):

#               Rc[iP, it] = -1.*((sigma_p*abs(SCC.fB(p, T))*xiGL)/fBA)
#               Rc_sh[iP, it] = np.nan
#               print("\n t < t_AB, sigmaAB/fBA looks like : ", Rc[iP, it])
          
#             elif t > SCC.tAB_RWSco(p):

#               Rc_sh[iP, it] = -1.*((sigma_p*abs(SCC.fB(p, T))*xiGL)/(-fBA))
#               Rc[iP, it] = np.nan        
          


# print(" \n Rc array looks like ",Rc, "\n Rc_sh array looks like ", Rc_sh)


#############################################################################
#####                rc loops with texture, evaluation                  #####
#############################################################################

for iP in range(0, len(Pressure), 1):
    print('\n\n Now P is:', Pressure[iP], '\n\n')

    p = Pressure[iP]

    sigma_p = tensionAB_p[iP]

    K_g = SCC.K(p) 
    
    for it in range(0, len(Temperature), 1):

       T = Temperature[it]

       t = T/SCC.Tcp(p)

       print('\n now temperature is:, ', T)

       if (p < 21.22) and (t <= 1.):

        GapA2 = (SCC.GapA(p, T))**2

        GapA = SCC.GapA(p, T)       

        fBA = SCC.fB(p, T) - SCC.fA(p, T) 

        # dimensional sigma_AB in SI unit
        sigma_p_fBxiGL = sigma_p*abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

        # thickness of A-phase texture, in unit of \xiGL(p, T), \xiGL(p)
        R = R_ratio*(SCC.xiGL_OC(p, T))
        # R = R_ratio*xiGL(p)
 
        Nomi = 4.*(R*fBA + sigma_p_fBxiGL)**2

        Dino1 = 135.*K_g*pi*pi*R*GapA2*(fBA*fBA)

        Dino2 = 8.*(R*fBA+sigma_p_fBxiGL)**(3.)

        Dino3 = 3.*np.sqrt(15.)*pi*GapA*fBA*np.sqrt(K_g*R*(Dino1 + 2*Dino2))

        Dino = np.cbrt(-Dino1-Dino2-Dino3)

        Nomi_Dino = Nomi/Dino

        fx = 4.*R*fBA-2.*sigma_p_fBxiGL+Nomi_Dino+Dino

        print("\n now p is ",p," fBA is ",fBA," Nomi is ",Nomi, " Dino1 is ", Dino1, " Dino2 is",Dino2," Dino3 is ",Dino3," Dino is ",Dino," Nomi/Dino is ",Nomi_Dino," fx is ",fx)
        
        # xiGL = SCC.xiGL_OC(p, T)

        fBA6_arr[iP, it] = 6.*fBA
        absfBA_xiGL_arr[iP, it] = abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

        R_term_arr[iP, it] = 4.*R*fBA
        sigma_term_arr[iP, it] = -2.*sigma_p_fBxiGL

        Nomi_Dino_arr[iP, it] = Nomi_Dino
        Dino_arr[iP, it] = Dino
        
        RcT[iP, it] = fx/(6.*fBA)
        RcT_sh[iP, it] = np.nan

        print(" \n pressure lower than P_PCP, r_c-Texture looks like, ", RcT[iP, it]*(10**(6.)))

       else: 
    
         if (t >= 1.):

           print(" bro, we just got temperature at Tc, save a np.nan. ")

           # mask data from t > 1 region

           fBA6_arr[iP, it] = np.nan
           absfBA_xiGL_arr[iP, it] = np.nan

           R_term_arr[iP, it] = np.nan
           sigma_term_arr[iP, it] = np.nan
           
           Nomi_Dino_arr[iP, it] = np.nan
           Dino_arr[iP, it] = np.nan

           RcT[iP, it] = np.nan

           RcT_sh[iP, it] = np.nan

      # only for A to B  
      # elif (t <1.) and (t >= h.t_AB(p)):

      # A to B or B to A
         elif (t == SCC.tAB_RWSco(p)):

            print(" bro, you got on T_AB, radius divergence. ")

            fBA6_arr[iP, it] = 6.*fBA
            absfBA_xiGL_arr[iP, it] = abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

            R_term_arr[iP, it] = 4.*R*fBA
            sigma_term_arr[iP, it] = -2.*sigma_p_fBxiGL


            Nomi_Dino_arr[iP, it] = np.nan
            Dino_arr[iP, it] = np.nan

            RcT[iP, it] = np.nan
            RcT_sh[iP, it] = np.nan
          
                 
         else:

            GapA2 = (SCC.GapA(p, T))**2.

            GapA = SCC.GapA(p, T)       

            fBA = SCC.fB(p, T) - SCC.fA(p, T) 
           
            # dimensional sigma_AB in SI unit
            sigma_p_fBxiGL = sigma_p*abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

            # thickness of A-phase texture, in unit of \xiGL(p, T), \xiGL(p)
            R = R_ratio*(SCC.xiGL_OC(p, T))
            # R = R_ratio*xiGL(p)
         
            Nomi = 4.*(R*fBA + sigma_p_fBxiGL)**2

            Dino1 = 135.*K_g*pi*pi*R*GapA2*(fBA*fBA)

            Dino2 = 8.*(R*fBA+sigma_p_fBxiGL)**(3.)

            Dino3 = 3.*np.sqrt(15.)*pi*GapA*fBA*np.sqrt(K_g*R*(Dino1 + 2*Dino2))

            Dino = np.cbrt(-Dino1-Dino2-Dino3)

            Nomi_Dino = Nomi/Dino

            fx = 4*R*fBA-2.*sigma_p_fBxiGL+Nomi_Dino+Dino

            print("\n now p is ",p," fBA is ",fBA," Nomi is ",Nomi, " Dino1 is ", Dino1, " Dino2 is",Dino2," Dino3 is ",Dino3," Dino is ",Dino," Nomi/Dino is ",Nomi_Dino," fx is ",fx)
        

            if t < SCC.tAB_RWSco(p):

              fBA6_arr[iP, it] = 6.*fBA
              absfBA_xiGL_arr[iP, it] = abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

              R_term_arr[iP, it] = 4.*R*fBA
              sigma_term_arr[iP, it] = -2.*sigma_p_fBxiGL

              Nomi_Dino_arr[iP, it] = Nomi_Dino
              Dino_arr[iP, it] = Dino

              RcT[iP, it] = fx/(6.*fBA)
              RcT_sh[iP, it] = np.nan
              print("\n t < t_AB, fx/(6.*fBA)) looks like : ", RcT[iP, it])
          
            elif t > SCC.tAB_RWSco(p):

              fBA6_arr[iP, it] = 6.*fBA
              absfBA_xiGL_arr[iP, it] = abs(SCC.fB(p, T))*(SCC.xiGL_OC(p, T))

              R_term_arr[iP, it] = 4.*R*fBA
              sigma_term_arr[iP, it] = -2.*sigma_p_fBxiGL

              Nomi_Dino_arr[iP, it] = Nomi_Dino
              Dino_arr[iP, it] = Dino

              RcT_sh[iP, it] = fx/(6.*fBA)
              RcT[iP, it] = np.nan        

# make the terms of r_c              
R_term_6fBA_arr = R_term_arr/fBA6_arr # 1st term
sigma_term_6fBA_arr = sigma_term_arr/fBA6_arr # 2nd term

NomiDino_6fBA_arr = Nomi_Dino_arr/fBA6_arr # 3rd term
Dino_6fBA_arr = Dino_arr/fBA6_arr # forth term

print(" \n Rc-Texture array looks like ",RcT)
print(" \n R_term_6fBA_arr is\n ", R_term_6fBA_arr, " \n sigma_term_6fBA_arr is\n ",sigma_term_6fBA_arr," \n NomiDino_6fBA_arr is\n ",NomiDino_6fBA_arr, "\n Dino_6fBA_arr\n ",Dino_6fBA_arr)



#############################################################################
#####      scatter plot the ratio verse the pressure                    #####
#############################################################################.

# fig, ax = plt.subplots(1,1)
# Levels1 = [0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2., 2.5, 3., 4., 5., 6.]
# Levels2 = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 2., 2.5, 3, 4., 5., 6., 7., 8., 9., 10., 11]

# Levels1sh = [0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9]

# # plot the TAB line with RWS2019 SC
# # TABrws = ax.contour(TAB_RWS.X*(10**3), TAB_RWS.Y, TAB_RWS.EnergyDensity_Difference_fABGL, levels=[0.0], colors='red')

# cf1 = ax.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc*(10**(6)), cmap=cm.PuBu_r, levels=Levels2)
# ax.set_ylabel(r'$p/bar$', fontsize=15.0,);
# ax.set_xlabel(r'$T$/mK', fontsize=15.0,);
# plt.xticks(fontsize=15.)
# plt.yticks(fontsize=15.)

# fig.colorbar(cf1, ax=ax, location = 'left')

# # cf2 = ax.contourf(h.T_GtoPlts6_poly(T_array*(10**(3)) + T_shift_array), P_array, Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels2)
# cf2 = ax.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels1sh)

# fig.colorbar(cf2, ax=ax, location = 'right')

# c3 = ax.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc*(10**(6)), levels=Levels1, colors='orange')
# plt.clabel(c3, inline=True, fontsize=10.5, colors='r')

# # c4 = ax.contour(h.T_GtoPlts6_poly(T_array*(10**(3)) +  T_shift_array), P_array, Rc_sh*(10**(6)), levels=Levels1, colors='blue')
# c4 = ax.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc_sh*(10**(6)), levels=Levels1sh, colors='blue')
# plt.clabel(c4, inline=True, fontsize=10.5, colors='b')

# # plot Greywall T_AB line
# # TABGreywall, = ax.plot(h.T_AB_Greywall_poly(Pressure), Pressure, color = "orange")

# # plot PLTS T_AB line from Greywall data
# # TABGreywalltoPlts, = ax.plot(h.T_GtoPlts6_poly(h.T_AB_Greywall_poly(Pressure)), Pressure, color = "red")

# # plot PLTS T_AB lien
# TABPlts, = ax.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "pink")

# # plot PLTS Tc
# # TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
# TCPlts, = ax.plot(h.T_GtoPlts6_low_poly(h.Tc_mK(Pressure)), Pressure, color = "purple")

# mSize = 140

# # scatter plot of parpia's constant pressure data
# sIC1 = ax.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label="p-T IC")
# sHEC1 = ax.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="p-T HEC")

# # scatter plot of parpia's constant Q data

# sHEC_CQ1 = ax.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="p-T HEC-CQ")
# sIC_CQ1 = ax.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="p-T IC-CQ")

# # ax.legend(labels=(r"$T_{AB}$", "p-T IC", "p-T HEC"))

# # leg1 = ax.legend(loc="upper right")
# # leg2 = ax.legend([TCPlts,TABGreywalltoPlts,sIC,sHEC],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC"],  fontsize=20.5, loc='lower right')
# leg2 = ax.legend([TCPlts,TABPlts,sIC1,sHEC1,sHEC_CQ1,sIC_CQ1],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC", r"$p-T-IC, ConsQ$",r"$p-T-HEC, ConsQ$"],  fontsize=15.0, loc='lower right')
# # ax.add_artist(leg1)

# # ax.set_title(r"contour plot of $(\sigma_{AB}(p)/f_{AB}(p,T))/{\mu}m$")
# ax.set_title(r"contour plot of $(-\sigma_{AB}(p)/f_{BA}(p,T))/{\mu}m$", fontsize=14.0,)
# ax.set_xlim([1.8, 2.4])
# ax.set_ylim([20., 30.])
# ax.grid(True)

# # fig.savefig("Rc_contour_RWS_GreywelltoPltsTC_AToBBToA5.pdf")

# plt.show()

#############################################################################
#####                contour plot of Nomi/Dino and Dino                 #####
#############################################################################

fig0, ax0 = plt.subplots(1,1)

c01 = ax0.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Nomi_Dino_arr/absfBA_xiGL_arr, cmap=cm.PuBu_r)
ax0.set_ylabel(r'$p/bar$', fontsize=15.0,);
ax0.set_xlabel(r'$T$/mK', fontsize=15.0,);
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)

fig0.colorbar(c01, ax=ax0, location = 'left')
ax0.set_title(r"$(Nomi/Dino)/|f_{BA}|\xi_{GL}$", fontsize=14.0)


fig2, ax2 = plt.subplots(1,1)

c21 = ax2.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Dino_arr/absfBA_xiGL_arr, cmap=cm.PuBu_r)
ax2.set_ylabel(r'$p/bar$', fontsize=15.0,);
ax2.set_xlabel(r'$T$/mK', fontsize=15.0,);
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)

fig2.colorbar(c21, ax=ax2, location = 'left')
ax2.set_title(r"$Dino/|f_{BA}|\xi_{GL}$", fontsize=14.0)


#############################################################################
#####                     plot every terms of r_c                     #######
#############################################################################

fig3, ax3 = plt.subplots(2,2)

# p = 23 bar
l301 = ax3[0,0].plot(Temperature*(10**3), R_term_6fBA_arr[6,:]*(10**6), color = "red", label=r"$\frac{2}{3}R$")
l302 = ax3[0,0].plot(Temperature*(10**3), sigma_term_6fBA_arr[6,:]*(10**6), color = "blue", label=r"$-\frac{\sigma_{AB}}{3{\Delta}f_{BA}}$")
l303 = ax3[0,0].plot(Temperature*(10**3), NomiDino_6fBA_arr[6,:]*(10**6), color = "cyan", label=r"$\frac{4(R{\Delta}f_{BA}+\sigma_{AB})^{2}}{6{\Delta}f_{BA}Dino}$")
l304 = ax3[0,0].plot(Temperature*(10**3), Dino_6fBA_arr[6,:]*(10**6), color = "green", label=r"$\frac{Dino}{6{\Delta}f_{BA}}$")

ax3[0,0].set_title(r"$p = 23 bar$")
ax3[0,0].set_ylim([-1, 5]) # micron
ax3[0,0].set_xlabel(r'$T/mK$', fontsize=15.0,)
ax3[0,0].set_ylabel(r'${\mu}m$', fontsize=15.0)

# leg00 = ax3[0,0].legend([l301,l302,l303,l304],[r"$\frac{2}{3}R$", r"$-\frac{\sigma_{AB}}{{\Delta}f_{BA}}$", r"3rd term of r_c", r"4th term of r_c"],  fontsize=15.0, loc='lower right')
leg00 = ax3[0,0].legend(loc="upper left", fontsize=15.0)
ax3[0,0].grid(True)


# p = 26 bar
l311 = ax3[0,1].plot(Temperature*(10**3), R_term_6fBA_arr[12,:]*(10**6), color = "red", label=r"$\frac{2}{3}R$")
l312 = ax3[0,1].plot(Temperature*(10**3), sigma_term_6fBA_arr[12,:]*(10**6), color = "blue", label=r"$-\frac{\sigma_{AB}}{3{\Delta}f_{BA}}$")
l313 = ax3[0,1].plot(Temperature*(10**3), NomiDino_6fBA_arr[12,:]*(10**6), color = "cyan", label=r"$\frac{4(R{\Delta}f_{BA}+\sigma_{AB})^{2}}{6{\Delta}f_{BA}Dino}$")
l314 = ax3[0,1].plot(Temperature*(10**3), Dino_6fBA_arr[12,:]*(10**6), color = "green", label=r"$\frac{Dino}{6{\Delta}f_{BA}}$")

ax3[0,1].set_title(r"$p = 26 bar$")
ax3[0,1].set_ylim([-1, 5]) # micron
ax3[0,1].set_xlabel(r'$T/mK$', fontsize=15.0)
ax3[0,1].set_ylabel(r'${\mu}m$', fontsize=15.0)

# leg01 = ax3[0,1].legend([l311,l312,l313,l314],[r"$\frac{2}{3}R$", r"$-\frac{\sigma_{AB}}{{\Delta}f_{BA}}$", r"3rd term of r_c", r"4th term of r_c"],  fontsize=15.0, loc='lower right')
leg01 = ax3[0,1].legend(loc="upper left", fontsize=15.0)
ax3[0,1].grid(True)


# p = 30 bar
l321 = ax3[1,0].plot(Temperature*(10**3), R_term_6fBA_arr[20,:]*(10**6), color = "red", label=r"$\frac{2}{3}R$")
l322 = ax3[1,0].plot(Temperature*(10**3), sigma_term_6fBA_arr[20,:]*(10**6), color = "blue", label=r"$-\frac{\sigma_{AB}}{3{\Delta}f_{BA}}$")
l323 = ax3[1,0].plot(Temperature*(10**3), NomiDino_6fBA_arr[20,:]*(10**6), color = "cyan", label=r"$\frac{4(R{\Delta}f_{BA}+\sigma_{AB})^{2}}{6{\Delta}f_{BA}Dino}$")
l324 = ax3[1,0].plot(Temperature*(10**3), Dino_6fBA_arr[20,:]*(10**6), color = "green", label=r"$\frac{Dino}{6{\Delta}f_{BA}}$")

ax3[1,0].set_title(r"$p = 30 bar$")
ax3[1,0].set_ylim([-1, 5]) # micron
ax3[1,0].set_xlabel(r'$T/mK$', fontsize=15.0)
ax3[1,0].set_ylabel(r'${\mu}m$', fontsize=15.0)

# leg10 = ax3[1,0].legend([l321,l322,l323,l324],[r"$\frac{2}{3}R$", r"$-\frac{\sigma_{AB}}{{\Delta}f_{BA}}$", r"3rd term of r_c", r"4th term of r_c"],  fontsize=15.0, loc='lower right')
leg10 = ax3[1,0].legend(loc="upper left", fontsize=15.0)
ax3[1,0].grid(True)

# p = 21 bar
l331 = ax3[1,1].plot(Temperature*(10**3), R_term_6fBA_arr[2,:]*(10**6), color = "red", label=r"$\frac{2}{3}R$")
l332 = ax3[1,1].plot(Temperature*(10**3), sigma_term_6fBA_arr[2,:]*(10**6), color = "blue", label=r"$-\frac{\sigma_{AB}}{3{\Delta}f_{BA}}$")
l333 = ax3[1,1].plot(Temperature*(10**3), NomiDino_6fBA_arr[2,:]*(10**6), color = "cyan", label=r"$\frac{4(R{\Delta}f_{BA}+\sigma_{AB})^{2}}{6{\Delta}f_{BA}(Dino}}$")
l334 = ax3[1,1].plot(Temperature*(10**3), Dino_6fBA_arr[2,:]*(10**6), color = "green", label=r"$\frac{Dino}{6{\Delta}f_{BA}}$")

ax3[1,1].set_title(r"$p = 21 bar$")
ax3[1,1].set_ylim([-1, 5]) # micron
ax3[1,1].set_xlabel(r'$T/mK$', fontsize=15.0)
ax3[1,1].set_ylabel(r'${\mu}m$', fontsize=15.0)

# leg11 = ax3[1,1].legend([l331,l332,l333,l334],[r"$\frac{2}{3}R$", r"$-\frac{\sigma_{AB}}{{\Delta}f_{BA}}$", r"3rd term of r_c", r"4th term of r_c"],  fontsize=15.0, loc='lower right')
leg11 = ax3[1,1].legend(loc="upper left", fontsize=15.0)
ax3[1,1].grid(True)




#############################################################################
#####    scatter plot the ratio verse the pressure  for Rc-Texture   ########
#############################################################################.

fig1, ax1 = plt.subplots(1,1)
Levels1 = [0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2., 2.5, 3., 4., 5., 6.]
Levels2 = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 2., 2.5, 3, 4., 5., 6., 7., 8., 9., 10., 11]

Levels1sh = [0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9]

# plot the TAB line with RWS2019 SC
# TABrws = ax.contour(TAB_RWS.X*(10**3), TAB_RWS.Y, TAB_RWS.EnergyDensity_Difference_fABGL, levels=[0.0], colors='red')

cf1 = ax1.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, RcT*(10**(6)), cmap=cm.PuBu_r, levels=Levels2)
ax1.set_ylabel(r'$p/bar$', fontsize=15.0,);
ax1.set_xlabel(r'$T$/mK', fontsize=15.0,);
plt.xticks(fontsize=15.)
plt.yticks(fontsize=15.)

fig1.colorbar(cf1, ax=ax1, location = 'left')

# cf2 = ax.contourf(h.T_GtoPlts6_poly(T_array*(10**(3)) + T_shift_array), P_array, Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels2)
#cf2 = ax.contourf(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc_sh*(10**(6)), cmap=cm.gist_heat, levels=Levels1sh)

#fig.colorbar(cf2, ax=ax, location = 'right')

c3 = ax1.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, RcT*(10**(6)), levels=Levels1, colors='orange')
plt.clabel(c3, inline=True, fontsize=10.5, colors='r')

# c4 = ax.contour(h.T_GtoPlts6_poly(T_array*(10**(3)) +  T_shift_array), P_array, Rc_sh*(10**(6)), levels=Levels1, colors='blue')
# c4 = ax.contour(h.T_GtoPlts6_low_poly(T_array*(10**(3))), P_array, Rc_sh*(10**(6)), levels=Levels1sh, colors='blue')
# plt.clabel(c4, inline=True, fontsize=10.5, colors='b')

# plot Greywall T_AB line
# TABGreywall, = ax.plot(h.T_AB_Greywall_poly(Pressure), Pressure, color = "orange")

# plot PLTS T_AB line from Greywall data
# TABGreywalltoPlts, = ax.plot(h.T_GtoPlts6_poly(h.T_AB_Greywall_poly(Pressure)), Pressure, color = "red")

# plot PLTS T_AB lien
TABPlts, = ax1.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "pink")

# plot PLTS Tc
# TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
TCPlts, = ax1.plot(h.T_GtoPlts6_low_poly(h.Tc_mK(Pressure)), Pressure, color = "purple")

mSize = 140

# scatter plot of parpia's constant pressure data
sIC1 = ax1.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label="p-T IC")
sHEC1 = ax1.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="p-T HEC")

# scatter plot of parpia's constant Q data

sHEC_CQ1 = ax1.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="p-T HEC-CQ")
sIC_CQ1 = ax1.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="p-T IC-CQ")

# ax.legend(labels=(r"$T_{AB}$", "p-T IC", "p-T HEC"))

# leg1 = ax.legend(loc="upper right")
# leg2 = ax.legend([TCPlts,TABGreywalltoPlts,sIC,sHEC],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC"],  fontsize=20.5, loc='lower right')
leg2 = ax1.legend([TCPlts,TABPlts,sIC1,sHEC1,sHEC_CQ1,sIC_CQ1],[r"$T_{c}^{PLTS}$", r"$T_{AB}^{PLTS}$", "p-T IC", "p-T HEC", r"$p-T-IC, ConsQ$",r"$p-T-HEC, ConsQ$"],  fontsize=15.0, loc='lower right')
# ax.add_artist(leg1)

# ax.set_title(r"contour plot of $(\sigma_{AB}(p)/f_{AB}(p,T))/{\mu}m$")
ax1.set_title(r"contour plot of $(fx(p,T)/6f_{BA}(p,T))/{\mu}m$", fontsize=14.0)
ax1.set_xlim([1.8, 2.4])
ax1.set_ylim([20., 30.])
ax1.grid(True)

# fig.savefig("Rc_contour_RWS_GreywelltoPltsTC_AToBBToA5.pdf")

plt.show()

