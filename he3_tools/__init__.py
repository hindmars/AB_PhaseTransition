#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 15:43:31 2022

he3_tools module

Functions for computing bulk free energy and other thermodynamic properties of 
equilibrium phases of superfluid He3.

@author: hindmars
"""

from .he3_constants import *
from .he3_data import *
from .he3_props import *
from .he3_matrix import *
from .he3_bases import *
from .he3_free_energy import *
from .he3_magnetic import *

for x in ["DEFAULT_SC_ADJUST", "DEFAULT_SC_CORRS", "DEFAULT_T_SCALE", "DEFAULT_ALPHA_TYPE"]:
    report_setting(x)


