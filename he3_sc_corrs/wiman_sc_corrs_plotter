#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:19:31 2022

@author: hindmars
"""

import matplotlib.pyplot as plt
import wiman_csvreader as w

plt.figure(figsize=(6,3))
for col in range(0,10,2):
    n = col//2+1
    plt.plot(w.data[:,col], w.data[:, col+1], marker='o', label=r'$\Delta\beta_{:d}$'.format(n))


plt.ylim(-0.5,0)
plt.xlim(0, 35)
plt.xlabel(r'$p$/bar')
plt.ylabel(r'$\Delta\beta^{\rm sc}_i/|\beta_1^{\rm wc}|$')
plt.grid()
plt.legend()
