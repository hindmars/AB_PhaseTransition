
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
#####                         plot contours                          ######
###########################################################################
    
# Levels1 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
Levels1 = np.arange(0.005, 1.02, 0.005) # micro meter
            
fig, ax = plt.subplots(1,2)

cf1 = ax[0].contourf(T_arr*(10**3), P_arr, xiGL_OC_arr*(10**6), cmap=cm.PuBu_r, levels=Levels1)

# fig.colorbar(cf1, ax=ax, location = 'left')

c1 = ax[0].contour(T_arr*(10**3), P_arr, xiGL_OC_arr*(10**6), levels=Levels1, colors='pink')
plt.clabel(c1, inline=True, fontsize=11, colors='r')

# plot Greywall T_AB line in axx[0]
TABGreywall1, = ax[0].plot(h.TAB_poly_Greywall(Pressure-h.p_pcp_bar), Pressure, color = "orange")
TcGreywall1, = ax[0].plot(h.Tc_poly_Greywall(Pressure), Pressure, color = "red")

leg1 = ax[0].legend([TABGreywall1,TcGreywall1],[r"$T_{AB}^{Greywall}$",r"$T_{c}^{Greywall}$"],  fontsize=15.0, loc='lower right')

ax[0].set_ylabel(r'$p/bar$', fontsize=15.0);
ax[0].set_xlabel(r'$T$/mK', fontsize=15.0);
ax[0].set_title(r"$\xi_{GL}^{OC}(p, T)/{\mu}m$, around PCP", fontsize=14.0)
ax[0].grid(True)

###########################################################################

cf2 = ax[1].contourf(T_arr*(10**3), P_arr, xiGL_RWS_arr*(10**6), cmap=cm.PuBu_r, levels=Levels1)

fig.colorbar(cf2, ax=ax, location = 'right')

c2 = ax[1].contour(T_arr*(10**3), P_arr, xiGL_RWS_arr*(10**6), levels=Levels1, colors='pink')
plt.clabel(c2, inline=True, fontsize=11, colors='r')


# plot Greywall T_AB line in axx[0]
TABGreywall2, = ax[1].plot(h.TAB_poly_Greywall(Pressure-h.p_pcp_bar), Pressure, color = "orange")
TcGreywall2, = ax[1].plot(h.Tc_poly_Greywall(Pressure), Pressure, color = "red")

# scatter plot of parpia's constant pressure data
sIC = ax[1].scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", label="p-T IC")
sHEC = ax[1].scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", label="p-T HEC")

leg2 = ax[1].legend([TABGreywall1,TcGreywall1, sIC, sHEC],[r"$T_{AB}^{Grewywall}$",r"$T_{c}^{Greywall}$",r"$p-T IC$",r"$p-T HEC$"],  fontsize=15.0, loc='lower right')

ax[1].set_ylabel(r'$p/bar$', fontsize=15.0);
ax[1].set_xlabel(r'$T$/mK', fontsize=15.0);
ax[1].set_title(r"$\xi_{GL}^{RWS19}(p, T)/{\mu}m$, around PCP", fontsize=14.0)
ax[1].grid(True)


plt.show()
