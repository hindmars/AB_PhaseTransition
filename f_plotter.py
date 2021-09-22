#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 13:35:33 2021

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h


mint, maxt = (0.7, 1.1)

t = np.linspace(mint, maxt, 100)


for p in range(0, 36, 4):

    plt.plot(t, h.f_A_norm(t, p) - h.f_B_norm(t, p), color='b', alpha=(p+4)/(34+4), label='{}'.format(p) )


plt.xlabel(r'$T/T_c$')
plt.ylabel(r'$\Delta f_{AB} /N(0) (k_B T_c)^2$')


plt.grid()
plt.xlim(mint,maxt+0.1)
plt.ylim(-1e-7, 1e-7)
plt.legend(title=r'$p$/bar')