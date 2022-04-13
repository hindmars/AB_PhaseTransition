
'''
contour plot of p-T GL coherent length

author: Quang, (timohyva@github)
'''


import Module_SC_Beta_V04 as SCC

import he3_tools_Vn01 as h


import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

import Module_New_Parpia_pT_parameter as Parpia
import Module_Parpia_pT_ConsQ as ConsQ



Pressure = np.arange(18.2,30.2,0.1) # bar
# Pressure = np.arange(20.0,24.01,0.01) # bar

Temperature = (np.arange(1.5, 2.5, 0.01))*(10**(-3)) # Kelvin
# Temperature = (np.arange(2.0, 2.4, 0.005))*(10**(-3)) # Kelvin

T_arr, P_arr = np.meshgrid(Temperature, Pressure)
print("\n T_arr looks like \n", T_arr,"\n P_arr looks like \n", P_arr )

xiGL_OC_arr = np.zeros(T_arr.shape)

xiGL_RWS_arr = np.zeros(T_arr.shape)



for ip in np.arange(0, len(Pressure), 1):

    p = Pressure[ip]
    
    for iT in np.arange(0, len(Temperature), 1):

        T = Temperature[iT]

        t = T/SCC.Tcp(p)

        print("\n now p is ", p," T is ", T, " t is ",t, "Tcp is ", SCC.Tcp(p))
        
        if t >= 1.:

            print("\n bro, T > Tc now, save a np.nan")

            xiGL_OC_arr[ip, iT] = np.nan

            xiGL_RWS_arr[ip, iT] = np.nan
            
        else:
        
            xiGL_OC_arr[ip, iT] = SCC.xiGL_OC(p, T)

            xiGL_RWS_arr[ip, iT] = SCC.xiGL_JWS(p, T)

    print("\n now xi_OC_arr looks like \n", xiGL_OC_arr*(10**6))

    print("\n now xi_RWS_arr looks like \n", xiGL_RWS_arr*(10**6))

###########################################################################
####                         mask TAB tails                          ######
###########################################################################

TAB_PLTS_masked = np.array([])

for p in Pressure:

    if p >= h.p_pcp_bar:
       TAB_PLTS_masked = np.append(TAB_PLTS_masked, h.TAB_poly_PLTS(p))

    elif p < h.p_pcp_bar:
       TAB_PLTS_masked = np.append(TAB_PLTS_masked, np.nan)        


    
###########################################################################
#####                         plot contours                          ######
###########################################################################
    
Levels1 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3]
# Levels1 = np.arange(0.005, 1.02, 0.005) # micro meter

mSize = 120 # scatter marker size
            
fig1, ax = plt.subplots(1,1)

cf1 = ax.contourf(h.T_GtoPlts6_low_poly(T_arr*(10**3)), P_arr, xiGL_OC_arr*(10**6), cmap=cm.PuBu_r, levels=Levels1)

# fig1.colorbar(cf1, ax=ax, location = 'left')

c1 = ax.contour(h.T_GtoPlts6_low_poly(T_arr*(10**3)), P_arr, xiGL_OC_arr*(10**6), levels=Levels1, colors='pink')
plt.clabel(c1, inline=True, fontsize=11, colors='r')

# # plot Greywall T_AB line in axx[0]
# TABGreywall1, = ax.plot(h.TAB_poly_Greywall(Pressure-h.p_pcp_bar), Pressure, color = "orange")
# TcGreywall1, = ax.plot(h.Tc_poly_Greywall(Pressure), Pressure, color = "red")

# plot PLTS T_AB lien
# TABPlts1, = ax.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "orange")
TABPlts1, = ax.plot(TAB_PLTS_masked, Pressure, color = "orange")

# plot PLTS Tc
# TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
TCPlts1, = ax.plot(h.T_GtoPlts6_low_poly(h.Tc_mK(Pressure)), Pressure, color = "red")


# scatter plot of parpia's constant pressure data
sIC1 = ax.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label="p-T IC")
sHEC1 = ax.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="p-T HEC")

# scatter plot of parpia's constant Q data

sHEC_CQ1 = ax.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="p-T HEC-CQ")
sIC_CQ1 = ax.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="p-T IC-CQ")


# leg1 = ax.legend([TABGreywall1,TcGreywall1, sIC1, sHEC1, sIC_CQ1, sHEC_CQ1],[r"$T_{AB}^{Greywall}$",r"$T_{c}^{Greywall}$",r"$p-T-IC$",r"$p-T-HEC$",r"$p-T-IC-CQ$",r"$p-T-HEC-CQ$"],  fontsize=15.0, loc='lower right')
leg1 = ax.legend([TABPlts1,TCPlts1, sIC1, sHEC1, sIC_CQ1, sHEC_CQ1],[r"$T_{AB}^{PLTS}$",r"$T_{c}^{PLTS}$",r"$p-T-IC$",r"$p-T-HEC$",r"$p-T-IC-CQ$",r"$p-T-HEC-CQ$"],  fontsize=15.0, loc='lower right')

ax.set_ylabel(r'$p/bar$', fontsize=15.0);
ax.set_xlabel(r'$T$/mK', fontsize=15.0);
ax.set_title(r"$\xi_{GL}^{OC}(p, T)/{\mu}m$, around PCP", fontsize=14.0)
ax.grid(True)

############################################################################

fig2, ax1 = plt.subplots(1,1)

cf2 = ax1.contourf(h.T_GtoPlts6_low_poly(T_arr*(10**3)), P_arr, xiGL_RWS_arr*(10**6), cmap=cm.PuBu_r, levels=Levels1)

# fig2.colorbar(cf2, ax=ax, location = 'left')

c2 = ax1.contour(h.T_GtoPlts6_low_poly(T_arr*(10**3)), P_arr, xiGL_RWS_arr*(10**6), levels=Levels1, colors='pink')
plt.clabel(c2, inline=True, fontsize=11, colors='r')


# # plot Greywall T_AB line in axx[0]
# TABGreywall2, = ax1.plot(h.TAB_poly_Greywall(Pressure-h.p_pcp_bar), Pressure, color = "orange")
# TcGreywall2, = ax1.plot(h.Tc_poly_Greywall(Pressure), Pressure, color = "red")

# plot PLTS T_AB lien
# TABPlts2, = ax1.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "orange")
TABPlts2, = ax1.plot(TAB_PLTS_masked, Pressure, color = "orange")

# plot PLTS Tc
# TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
TCPlts2, = ax1.plot(h.T_GtoPlts6_low_poly(h.Tc_mK(Pressure)), Pressure, color = "red")

# scatter plot of parpia's constant pressure data
sIC2 = ax1.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label="p-T IC")
sHEC2 = ax1.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="p-T HEC")

# scatter plot of parpia's constant Q data
sHEC_CQ2 = ax1.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="p-T HEC-CQ")
sIC_CQ2 = ax1.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="p-T IC-CQ")


# leg2 = ax1.legend([TABGreywall2,TcGreywall2, sIC2, sHEC2, sIC_CQ2, sHEC_CQ2],[r"$T_{AB}^{Grewywall}$",r"$T_{c}^{Greywall}$",r"$p-T-IC$",r"$p-T-HEC$",r"$p-T-IC-CQ$",r"$p-T-HEC-CQ$"],  fontsize=15.0, loc='lower right')
leg2 = ax1.legend([TABPlts2,TCPlts2, sIC2, sHEC2, sIC_CQ2, sHEC_CQ2],[r"$T_{AB}^{PLTS}$",r"$T_{c}^{PLTS}$",r"$p-T-IC$",r"$p-T-HEC$",r"$p-T-IC-CQ$",r"$p-T-HEC-CQ$"],  fontsize=15.0, loc='lower right')

ax1.set_ylabel(r'$p/bar$', fontsize=15.0);
ax1.set_xlabel(r'$T$/mK', fontsize=15.0);
ax1.set_title(r"$\xi_{GL}^{RWS19}(p, T)/{\mu}m$, around PCP", fontsize=14.0)
ax1.grid(True)


plt.show()
