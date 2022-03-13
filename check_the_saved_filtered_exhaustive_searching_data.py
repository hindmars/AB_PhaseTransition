

import numpy as np
import matplotlib.pyplot as plot1


with np.load('filtered_lambdaBar_Dalphai_8200loops.npz') as data:
    d11 = data['d11']
    d12 = data['d12']
    d13 = data['d13']
    d21 = data['d21']
    d22 = data['d22']
    d23 = data['d23']
    d31 = data['d31']
    d32 = data['d32']
    d33 = data['d33']

    psi11 = data['psi11']
    psi12 = data['psi12']
    psi13 = data['psi13']
    psi21 = data['psi21']
    psi22 = data['psi22']
    psi23 = data['psi23']
    psi31 = data['psi31']
    psi32 = data['psi32']
    psi33 = data['psi33']
    
    lambda_bar = data['lambdaBar']

    M2 = data['M2']
    delta = data['delta']
    lambda_pT = data['lambda_pT']
    
    print(d11, " \n and the type is :", d11.dtype, " and the size is : ", d11.size)
    print(d32, " \n and the type is :", d32.dtype, " and the size is : ", d32.size)
    print(lambda_bar, " \n and the type is :", lambda_bar.dtype, " and the size is : ", lambda_bar.size)


print(" \n\n show the nagative elements of M2 : \n ",M2[M2 < 0], " \n\n show the nagative elements of \delta : \n ",delta[delta < 0], " \n\n show the nagative elements of \lambda : \n ",lambda_pT[lambda_pT < 0])


print(" the number of nagative elements of \delta is : ", (delta[delta < 0]).size)

# positive lambda filter
positive_lambda_filter = delta > 0
print(" \n\n the positive lambda filter is : \n ", positive_lambda_filter)

lambda_bar_filtered = lambda_bar[positive_lambda_filter]

d11_filtered = d11[positive_lambda_filter]
d12_filtered = d12[positive_lambda_filter]
d13_filtered = d13[positive_lambda_filter]
d21_filtered = d21[positive_lambda_filter]
d22_filtered = d22[positive_lambda_filter]
d23_filtered = d23[positive_lambda_filter]
d31_filtered = d31[positive_lambda_filter]
d32_filtered = d32[positive_lambda_filter]
d33_filtered = d33[positive_lambda_filter]


psi11_filtered = psi11[positive_lambda_filter]
psi12_filtered = psi12[positive_lambda_filter]
psi13_filtered = psi13[positive_lambda_filter]
psi21_filtered = psi21[positive_lambda_filter]
psi22_filtered = psi22[positive_lambda_filter]
psi23_filtered = psi23[positive_lambda_filter]
psi31_filtered = psi31[positive_lambda_filter]
psi32_filtered = psi32[positive_lambda_filter]
psi33_filtered = psi33[positive_lambda_filter]

# polar plot these D_alphai with positive \delta
markersaize = 3.6
fig1, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plot1.subplots(3, 3, subplot_kw=dict(projection='polar'))

ax1.scatter(psi11_filtered, d11_filtered, s = markersaize)
ax2.scatter(psi12_filtered, d12_filtered, s = markersaize)
ax3.scatter(psi13_filtered, d13_filtered, s = markersaize)
ax4.scatter(psi21_filtered, d21_filtered, s = markersaize)
ax5.scatter(psi22_filtered, d22_filtered, s = markersaize)
ax6.scatter(psi23_filtered, d23_filtered, s = markersaize)
ax7.scatter(psi31_filtered, d31_filtered, s = markersaize)
ax8.scatter(psi32_filtered, d32_filtered, s = markersaize)
ax9.scatter(psi33_filtered, d33_filtered, s = markersaize)

# plot the filtered lambda_bar
fig2, ax_lambda_bar = plot1.subplots(1, 1)
ax_lambda_bar.scatter(np.arange(0,lambda_bar_filtered.size,1), lambda_bar_filtered);ax_lambda_bar.grid()

plot1.show()
