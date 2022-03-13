












import csv

# import the csv data of experment
with open('slow_transition_lot_et_al_2021.csv', newline='') as f:
     reader = csv.reader(f)
     data = list(reader)

# print(data)


list1 = list(zip(*data))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

Temperature = []
for ii in list(list1[0]):
    Temperature.append(float(ii))

print("Lotynk's Temperature is", Temperature)  
    
pressure = []
for ii in list(list1[1]):
    pressure.append(float(ii))    

print("Lotynk's pressure is", pressure)  
