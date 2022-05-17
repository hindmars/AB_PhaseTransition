'''
This script is for trying out the good fudge coeffecient arr for beta_i

13th. April. 2022

author: Quang (timohyva@github)
'''


import numpy as np
import math

import Module_critical_radius_AB_wall_ContourPlot_V01 as crcAB
import Module_SC_Beta_V05 as SCC


# sample_number = 10
sample_number = 1
#rws_arr = np.array([1., 1.08, 1., 1., 1.])
# rws_arr = np.array([1., 1.18, 1., 0.86, 0.86])
# rws_arr = np.array([1.1, 1.20, 1.1, 0.86, 0.86])
rws_arr = np.array([1.8, 1.2, 1.8, 0.98, 0.98])

# dfudge_Arr = np.random.uniform(-0.1, 0.1, size=(sample_number, 5))
# df_Arr_245 = np.random.uniform(-0.1, 0.1, size=(sample_number, 3))
df_Arr_245 = np.random.uniform(-0.2, -0.1, size=(sample_number, 2))
# print("\n df_arr_245 looks like\n ", df_Arr_245)

a = np.zeros((sample_number, 1))
# dfudge_Arr = np.concatenate((a,df_Arr_245[:,0:1],a,df_Arr_245[:,1:2],df_Arr_245[:,2:3]),axis=1)
dfudge_Arr = np.concatenate((a,a,a,df_Arr_245[:,0:1],df_Arr_245[:,1:2]),axis=1)

print("\n defugde_Arr looks like\n ", dfudge_Arr)


for n in np.arange(0, sample_number, 1):

    print("\n now n is", n, ", dfudge_arr looks like\n ", dfudge_Arr[n,:])

    # fudge_arr = rws_arr + dfudge_Arr[n,:]
    fudge_arr = rws_arr 

    print("\n now n is", n, ", fudge_arr looks like\n ", fudge_arr, "\n\n")

    # SCC.fc_arr = fudge_arr

    # print("\n now n is ",n,", SCC.fc_arr is ", SCC.fc_arr, "\n\n\n" )

    crcAB.plot_contour(fudge_arr)
