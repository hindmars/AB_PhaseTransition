












import pandas as pd
import numpy as np

# import the csv data of Parpia's new data
sheet = pd.read_excel('Const_P.xlsx')
sheet = np.array(sheet)


# pressure
pressure = sheet[:, 1]; print("\n\n Parpia's new pressure is ",pressure)

# T_IC, unit mK
T_IC = sheet[:, 2]; print("\n Parpia's T_IC is ", T_IC)

# T_HEC
T_HEC = sheet[:, 3]; print("\n Parpia's T_HEC is ", T_HEC)

