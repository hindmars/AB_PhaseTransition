#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 08:48:20 2022

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import nucleation_tools as nt


fig, ax = plt.subplots(2,1, figsize=(5,5))

T = np.linspace(0.5,1,500)
Tn = 0.7
eta = 100


ax[0].plot(T, nt.nuc_prob_dist(T, Tn, eta) )

ax[0].set_ylabel(r'$dP/dT$')

ax[0].grid(True)
ax[0].set_title(r'Probability distribution: $T_f/T_{{AB}} = {:.2f}$, $|\eta_f| = {:.0f}$'.format(Tn, eta))
ax[0].set_xlabel(r'$T/T_{\rm AB}$')
ax[0].set_ylabel(r'$dP/dT$')

ax[1].plot(T, nt.nuc_prob(T, Tn, eta) )
ax[1].set_xlabel(r'$T/T_{\rm AB}$')
ax[1].set_ylabel(r'$P$')
ax[1].grid(True)
ax[1].set_title('Cumulative nucleation probability (cooling)')
ax[1].set_xlim(min(T),max(T))

fig.tight_layout()

fig.savefig('nuc_prob.pdf')


