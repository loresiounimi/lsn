import numpy as np
import matplotlib.pyplot as plt

#per tutti e 11 i rank creo i grafici
for i in range (11):
    #carico i dati
    data = np.loadtxt(f"OUTPUT/best_path_cartesian_coordinates_{i}.dat", comments="#")#f per leggere i come variabile del ciclo for
    x = data[:,1]
    y = data[:,2]
    #aggiungo la prima città all'ultimo posto per chiudere il percorso nel grafico
    x = np.append(x, x[0])
    y = np.append(y, y[0])

    #creo il grafico
    plt.figure(figsize=(8, 8))
    plt.plot(x, y, '-o', color='blue', markersize=5, label='Percorso TSP')
    #evidenzio la città di partenza
    plt.plot(x[0], y[0], 'ro', markersize=10, label='Prima città')
    plt.text(x[0]+0.10, y[0] - 0.10, 'Start', color='red', ha='center', fontsize=12)
    plt.title(f"Percorso TSP rank {i}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"grafici/TSP_{i}.png")

    plt.close()

    #carico i dati
    data1 = np.loadtxt(f"OUTPUT/best_path_{i}.dat", comments="#")
    x1 = data1[:,0]
    y1 = data1[:,1]

    #creo il grafico
    plt.figure(figsize=(8, 8))
    plt.plot(x1, y1, '-', color='blue')
    plt.title(fr"$L_1(x)$ best in funzione della generazione rank {i}")#f per leggere i dal ciclo for, r per Latex
    plt.xlabel("Generazione")
    plt.ylabel(r"$L_1(x)$ best")
    plt.grid(True)
    plt.savefig(f"grafici/Percorso_migliore_{i}.png")

    plt.close()

    #carico i dati
    data2 = np.loadtxt(f"OUTPUT/best_path_half_population_{i}.dat", comments="#")
    x2 = data2[:,0]
    y2 = data2[:,1]

    #creo il grafico
    plt.figure(figsize=(8, 8))
    plt.plot(x2, y2, '-', color='blue')
    plt.title(fr"$\langle L_1(x) \rangle$ in funzione della generazione rank {i}")
    plt.xlabel("Generazione")
    plt.ylabel(r"$\langle L_1(x) \rangle$")
    plt.grid(True)
    plt.savefig(f"grafici/Percorso_migliore_half_population_{i}.png")
    plt.close()
