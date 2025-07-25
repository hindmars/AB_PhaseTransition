#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 12:06:34 2024

test sympy soling cubic

@author: hindmars
"""

from sympy.solvers import solve
from sympy import Symbol
from sympy import symbols
x = Symbol('x')

a, b, c, d = symbols('alpha beta gamma delta')

foo = solve(a + 2*b*x + 3*c*x**2 + 4*d*x**3, x)

foo_num_list = solve(-1 + 2*x - 3*x**2 + 4*x**3, x)


for foo_num in foo_num_list:
    print(foo_num.evalf())
