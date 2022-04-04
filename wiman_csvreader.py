#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:46 2022

J Wiman thesis sc corrections reader

@author: hindmars
"""

import numpy as np

data = np.genfromtxt('wiman-thesis_sc_corrs_3_sorted.csv', delimiter=',', skip_header=2)

names = ['p_nodes', 'c1', 'c2', 'c3', 'c4', 'c5']

dtype=[('p', 'float'), ('val', 'float')]

# for n in names:
#     dtype.append((n, 'float'))

print(dtype)

p_nodes = data[:,0]
c1 = data[:,1]
c2 = data[:,3]
c3 = data[:,5]
c4 = data[:,7]
c5 = data[:,9]

def pretty_print_list(xlist, name):
    print(name + ' = [', end='')
    print(*xlist, sep=", ", end='')
    print(']')



for liszt, name in zip([p_nodes, c1, c2, c3, c4, c5], names):
    pretty_print_list(liszt, name)



