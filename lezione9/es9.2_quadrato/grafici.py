import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data = np.loadtxt("OUTPUT/best_path_cartesian_coordinates.dat", comments="#")
x = data[:,1]
y = data[:,2]
#aggiungo la prima città all'ultimo posto per chiudere il percorso nel grafico
x = np.append(x, x[0])
y = np.append(y, y[0])

#creo il grafico
plt.figure(figsize=(8, 8))
#plotto la soluzione del TSP
plt.plot(x, y, '-o', color='blue', markersize=5, label='Percorso TSP')
#evidenzio la città iniziale e finale del percorso
plt.plot(x[0], y[0], 'ro', markersize=10, label='Prima città')
plt.text(x[0]+0.10, y[0] - 0.10, 'Start', color='red', ha='center', fontsize=12)
plt.title("Percorso TSP")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.grid(True)
plt.savefig("grafici/TSP.png")

plt.clf()

#carico i dati
data1 = np.loadtxt("OUTPUT/best_path.dat", comments="#")
x1 = data1[:,0]
y1 = data1[:,1]

#creo il grafico
plt.figure(figsize=(8, 8))
plt.plot(x1, y1, '-', color='blue')
plt.title(r"$L_1(x)$ best in funzione della generazione")
plt.xlabel("Generazione")
plt.ylabel(r"$L_1(x)$ best")
plt.grid(True)
plt.savefig("grafici/Percorso_migliore.png")

plt.clf()

#carico i dati
data2 = np.loadtxt("OUTPUT/best_path_half_population.dat", comments="#")
x2 = data2[:,0]
y2 = data2[:,1]

#creo il grafico
plt.figure(figsize=(8, 8))
plt.plot(x2, y2, '-', color='blue')
plt.title(r"$\langle L_1(x) \rangle$ in funzione della generazione")
plt.xlabel("Generazione")
plt.ylabel(r"$\langle L_1(x) \rangle$")
plt.grid(True)
plt.savefig("grafici/Percorso_migliore_half_population.png")

