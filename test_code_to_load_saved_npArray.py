import numpy as np


with np.load('filtered_lambdaBar_Dalphai.npz') as data:
    x2 = data['d11']
    x1 = data['d32']
    y2 = data['lambdaBar']
    print(x1, " \n and the type is :", x1.dtype)
    print(x2, " \n and the type is :", x2.dtype)
    print(y2, " \n and the type is :", y2.dtype)

