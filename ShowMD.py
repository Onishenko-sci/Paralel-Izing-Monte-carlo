from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as figure
from tkinter import *
import csv
import sys
from scipy.optimize import curve_fit

# Реализация возможности выбора файла для отрисовки. По умолчанию '../render/md_render.txt'
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = './' + sys.argv[1]
    else:
        filename = '../render/md_render.txt'


# Считывание информации из файла.
i = 0
j = 0
head_readed = False
read_corelation_data = False
with open(filename, newline='\n') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:

        if (head_readed == False):
            bound_x = float(row[0])
            bound_y = float(row[1])
            Number_of_particles = int(row[2])
            radius = float(row[3])
            steps = int(row[4])
            delta_t = float(row[5])
            frame_rate = int(row[6])
            cor_points = int(row[7])
            frames = int(steps/frame_rate)
            frame_x = np.zeros(int(Number_of_particles*(frames+1))
                               ).reshape(Number_of_particles, int(frames+1))
            frame_y = np.zeros(int(Number_of_particles*(frames+1))
                               ).reshape(Number_of_particles, int(frames+1))

            current_frame = np.zeros(int(frames+1))
            temperature = np.zeros(int(frames))
            potential_energy = np.zeros(int(frames))
            kinetic_energy = np.zeros(int(frames))
            sqared_displacment = np.zeros(int(frames))
            corl_x = np.zeros(int(row[7]))
            corl = np.zeros(int(row[7]))
            head_readed = True
            continue

        if ((i == Number_of_particles) and (not read_corelation_data)):
            current_frame[j] = float(row[1])
            kinetic_energy[j] = float(row[2])
            potential_energy[j] = float(row[3])
            sqared_displacment[j] = float(row[4])
            temperature[j] = float(row[5])
            i = 0
            j = j+1
            continue

        if (float(row[0]) == 101):
            read_corelation_data = True
            i = 0
            continue

        if (read_corelation_data):
            corl_x[i] = float(row[0])
            corl[i] = float(row[1])
            i = i+1
            continue

        frame_x[i, j] = float(row[0])
        frame_y[i, j] = float(row[1])
        i = i+1

# Апроксимация квадратичного отклонения


def sqr(x, a):
    return a*x**2


X = np.arange(0, len(sqared_displacment))*frame_rate
X = X*delta_t*(10**9)
Y = sqared_displacment/(radius)**2

args, covar = curve_fit(sqr, X, Y)
Y_predicted = sqr(X, args[0])


# Вывод графиков наблюдаемых
plt.figure(figsize=(11, 5))
plt.subplot(1, 2, 1)
plt.title('Correlation function')
plt.ylabel('Correlation function')
plt.xlabel('Distance between particles, Angstrom')
plt.plot(corl_x*(10**10), corl, 'b.')
plt.plot(np.array([radius, radius])*(10**10),
         np.linspace(0, corl.max(), 2), 'g--')  # Вывод на график радиуса частицы в виде вертикальной линии

plt.subplot(1, 2, 2)
plt.title('Mean squared displacement')
plt.ylabel('Mean squared displacement, sqared atom radius')
plt.xlabel('Time, ns')
plt.plot(X, Y, "b.")
plt.plot(X, Y_predicted, "r--")

plt.show()

total_e = kinetic_energy+potential_energy

plt.figure(figsize=(10, 7), dpi=120)
plt.subplot(2, 1, 1)
plt.title('Temperature')
plt.ylabel('Temperature, K')
plt.plot(X, temperature, 'b.')
plt.subplot(2, 1, 2)
plt.title('Energy')
plt.ylabel('Energy, J')
plt.plot(X, potential_energy, 'b-', label="Potential")
plt.plot(X, kinetic_energy, 'r-', label='Kinetic')
plt.plot(X, total_e, 'g-', label='Total')
plt.xlabel('Time, ns')
plt.legend()
plt.show()

# Константы для отрисовки
mashtab = 800/bound_x
scene_x = 5
scene_y = 5

# Радиус частицы уменьшен в 2 раза для наглядности
radius = radius/2
density = (Number_of_particles*radius**2)/(bound_x*bound_y)

# Траектория отдельной частицы
traked = 0  # Номер отслеживаемой частицы
tragectory_x = frame_x[traked, :]*mashtab
tragectory_y = frame_y[traked, :]*mashtab

# Создание окна
root = Tk()
c = Canvas(root, width=int(bound_x*mashtab+500),
           heigh=int(bound_y*mashtab+30), bg="#303030", highlightthickness=0)
c.pack()

# Отрисовка траектории частицы
for k in range(int(frames)):
    c.create_oval(tragectory_x[k]-2+5, tragectory_y[k]-2+5,
                  tragectory_x[k]+2+5, tragectory_y[k]+2+5, width=0, fill="white")

# Отрисовка границ, неизменяемых переменных.
c.create_line(bound_x*mashtab+5, 5, bound_x*mashtab+5,
              bound_y*mashtab+5, width=5, fill="#BDD4F1")
c.create_line(bound_x*mashtab+5, bound_y*mashtab+5, 5,
              bound_y*mashtab+5, width=5, fill="#BDD4F1")
c.create_line(5, 5, 5, bound_y*mashtab+5, width=5, fill="#BDD4F1")
c.create_line(5, 5, bound_x*mashtab+5, 5, width=5, fill="#BDD4F1")
c.create_text(bound_x*mashtab+10, bound_y*mashtab-15, text=("Particles: " +
              str(Number_of_particles)), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'))
c.create_text(bound_x*mashtab+10, bound_y*mashtab-45, text=("Density: " +
              str(density)[0:5]), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'))

# Основной цикл отрисовки
for i in range(frames):
    c.delete("del")
    # Вывод наблюдаемых
    c.create_text(bound_x*mashtab+10, 5,   text=("Time: " + str(current_frame[i]) + '*' + str(
        delta_t)), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10, 35,  text=("Temperature: " + str(
        temperature[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10, 65,  text=("Kinetic Energy: " + str(
        kinetic_energy[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10, 95,  text=("Potential Energy: " + str(
        potential_energy[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10, 125,  text=("Total Energy: " + str(
        potential_energy[i]+kinetic_energy[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10, 155, text=("Mean squared displacement: " + str(
        sqared_displacment[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")

    # Вывод индикатора выполнения (линия внизу окна)
    c.create_line(0, bound_y*mashtab+25, (int(i/(frames/1000))/1000)*(bound_x *
                  mashtab+500), bound_y*mashtab+25, width=10, fill="#6BAABF", tags="del")

    # Вывод частиц
    for j in range(Number_of_particles):
        if j == traked:
            c.create_oval(scene_x + (frame_x[j, i]-radius)*mashtab,
                          scene_y + (frame_y[j, i]-radius)*mashtab,
                          scene_x + (frame_x[j, i]+radius)*mashtab,
                          scene_y + (frame_y[j, i]+radius)*mashtab,
                          fill='#FB0000', tags="del", activefill="blue", outline="#FB0000")
        else:
            c.create_oval(scene_x + (frame_x[j, i]-radius)*mashtab,
                          scene_y + (frame_y[j, i]-radius)*mashtab,
                          scene_x + (frame_x[j, i]+radius)*mashtab,
                          scene_y + (frame_y[j, i]+radius)*mashtab,
                          fill="#6BAABF", tags="del", outline="#6BAABF")

    root.update()

    # Консольный индикатор выполнения.
    if (i % (frames/10) == 0):
        print(int(i/(frames/100)), '%')
        if (int(i/(frames/100)) == 99):
            print("Done!")

root.mainloop()
