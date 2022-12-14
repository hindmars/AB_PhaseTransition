#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 10:03:06 2022

Fit Hakonen et al 1983 nucelation data to nucleation theory

@author: hindmars
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import nucleation_tools as nt

#%%

fig, ax = plt.subplots()

data = np.loadtxt('Hak+83_histo_outline.csv', delimiter=',')

t_bin_edges = data[:,0]
events_edges = data[:,1]

t_pre = np.arange(0.6, 0.655, 0.005)
t_post = np.arange(0.720, 0.755, 0.005)

events_edges_pre = np.zeros_like(t_pre)
events_edges_post = np.zeros_like(t_post)

t_bin_edges = np.concatenate((t_pre, t_bin_edges,t_post))
events_edges = np.concatenate((events_edges_pre, events_edges, events_edges_post))

#%%

t = (t_bin_edges[:-1] + t_bin_edges[1:])/2
dt = -np.mean((t_bin_edges[:-1] - t_bin_edges[1:]))
events = (events_edges[:-1] + events_edges[1:])/2

# plt.plot(t_bin_edges, events_edges)

# plt.plot(t, events, 'k.', label='Hakonen et al 1983')
plt.plot(t_bin_edges, events_edges, 'k', label='Hakonen et al 1983')

#**

popt, cov = opt.curve_fit(nt.nuc_prob_dist, t, events, p0=(0.67,100,0.01))

print(popt)

#%%

t_smooth = np.linspace(min(t), max(t),500)

plt. plot(t_smooth, nt.nuc_prob_dist(t_smooth, *popt), label=r'$dP/dT, T_f = {:.3f}$, $\eta_f = {:.1f}$'.format(*popt))

plt.xlabel(r'$T/T_{\rm c}$')
plt.ylabel(r'$N$')
plt.xlim(0.6,0.76)
plt.ylim(0,5)
plt.title(r'Fit to nucleation data with $P = \exp( - |T_f/T|^{\eta_f} )$')
plt.legend()

plt.savefig('Hak+83_fit.pdf')