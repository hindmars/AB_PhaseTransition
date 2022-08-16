












import csv

# import the csv data of experment
with open('parpia_fig3b_HEC.csv', newline='') as f1:
     reader1 = csv.reader(f1)
     data1 = list(reader1)

# import the csv data of experment
with open('parpia_fig3b_IC.csv', newline='') as f2:
     reader2 = csv.reader(f2)
     data2 = list(reader2)
     
# print(data)


list1 = list(zip(*data1))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

T_HEC = []
for ii in list(list1[0]):
    T_HEC.append(float(ii))

print("\nConsQ HEC's Temperature is\n", T_HEC)  
    
p_HEC = []
for ii in list(list1[1]):
    p_HEC.append(float(ii))    

print("\nConsQ HEC's pressure is\n", p_HEC)  


list2 = list(zip(*data2))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

T_IC = []
for ii in list(list2[0]):
    T_IC.append(float(ii))

print("\nConsQ IC's Temperature is\n", T_IC)  
    
p_IC = []
for ii in list(list2[1]):
    p_IC.append(float(ii))    

print("\nConsQ IC's pressure is\n", p_IC)  
