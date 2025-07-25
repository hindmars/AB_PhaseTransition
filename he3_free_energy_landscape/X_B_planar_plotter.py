#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:25:44 2022

@author: hindmars
"""


import numpy as np
import matplotlib.pyplot as plt
import he3_tools as h
from free_energy_section_plot import *


p = 25.

for t in [0, 0.4, 0.8]:
    plot_free_energy_sections(t, p, savefig=False)

