import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as figure


input_ar = np.zeros((20,2))

fr = np.array([0.5,1.0,2.5,5.0,10.0])
i=0
j=0
filename="fr_rate.txt"
with open(filename, newline='\n') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        input_ar[i,0] = float(row[0]) #solid block
        input_ar[i,1] = float(row[1]) #


        i+=1

plt.figure(figsize=(15, 5), dpi=120)
print(input_ar)

plt.title('Frame rate')
plt.ylabel('Time, sec')
plt.xlabel('Threads')
plt.ylim(0,6)
counter=0
for i in fr:
    plt.plot(np.array([1,2,3,4]), input_ar[counter,:] , 'o-', label=str(i))
    counter+=1
    
plt.legend()


plt.show()