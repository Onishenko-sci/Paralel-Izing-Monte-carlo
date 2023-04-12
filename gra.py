import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as figure


input_ar = np.zeros((6,4))


Temperature = np.array([0.5,1.0,2.5,5.0,10.0])
with_block = np.zeros((2,5))
sim = np.zeros((2,5))
fr = np.array([2,4,16,64,128, 0])
i=0
j=0
filename="fr_rate.txt"
with open(filename, newline='\n') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        with_block[0,i] = float(row[0])
        sim[0,i]  = float(row[1]) 
        with_block[1,i]  = float(row[2]) 
        sim[1,i] = float(row[3]) 

        i+=1

plt.figure(figsize=(15, 5), dpi=120)
print(input_ar)

#plt.title('Зависимость времени раcxtnf от частоты использования барьера')
plt.ylabel('Время расчета, сек')
plt.xlabel('Температура')
#plt.xticks([1,2,3,4])
#plt.ylim(0,5)
counter=0

#plt.plot(Temperature, with_block[1,:] , 'o-', label="block4")
plt.plot(Temperature, with_block[0,:] , 'o-', label="block2")
#plt.plot(Temperature, sim[1,:] , 'o-', label="sim4")
plt.plot(Temperature, sim[0,:] , 'o-', label="Раздельная блокировка 2")

    
plt.legend()


plt.show()