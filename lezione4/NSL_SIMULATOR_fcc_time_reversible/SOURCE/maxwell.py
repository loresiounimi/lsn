import numpy as np
import matplotlib.pyplot as plt

#nella simulazione backward sono stati presi i dati per più blocchi (a parità di step) rispetto alla simulazione forward, questo per verificare che effettivamente il sistema non torna mai alla configurazione iniziale di bassa entropia, anzi, dopo un iniziale riordinamento, il sistema tende a una configurazione di alta entropia. 
#fenomeno dovuto all'aritmetica e precisione finita della macchina, non a verlet. all'aumentare della simulazione i piccoli errori dovuti alla macchina si accumulano e provocano un allontanamento significativo e divergente dalla simulazione forward. 
#nella simulazione forward sono stati presi blocchi non fino all'equilibrio proprio per mettere in evidenza la non reversibilità di verlet che ci si aspetta teoricamente. è più evidente la parziale simmetria temporale se nella backward simulation non si parte da una situazione di equilibrio (come in questo caso)

# Carica tutti i dati
dataforwardT = np.loadtxt("../OUTPUT/FORWARD_TIME/temperature.dat", comments="#")
databackwardT = np.loadtxt("../OUTPUT/temperature.dat", comments="#")
dataforwardK = np.loadtxt("../OUTPUT/FORWARD_TIME/kinetic_energy.dat", comments="#")
databackwardK = np.loadtxt("../OUTPUT/kinetic_energy.dat", comments="#")
dataforwardP = np.loadtxt("../OUTPUT/FORWARD_TIME/potential_energy.dat", comments="#")
databackwardP = np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")

# Dati per il forward (asse x traslato di -7)
x_forwardT = dataforwardT[:, 0] - 7
y_forwardT = dataforwardT[:, 1]

# Dati per il backward (asse x traslato di -1)
x_backwardT = databackwardT[:, 0] - 1
y_backwardT = databackwardT[:, 1]

# Plot 
plt.plot(x_forwardT, y_forwardT, label="Forward simulation", color="blue", marker='o', markersize=3)
plt.plot(x_backwardT, y_backwardT, label="Backward simulation", color="orange", marker='o', markersize=3)
plt.xlabel("Direzione del tempo")
plt.ylabel("Temperatura")
plt.title('Verifica della Reversibilità Temporale: Andata e Ritorno della Temperatura' )
plt.grid(True)
plt.text(0.59, 0.7,"Nota: in 0 finisce la \nsimulazione forward e inizia\nla simulazione backward", transform=plt.gca().transAxes, fontsize=9, fontweight='bold')#aggiungo un commento al grafico, non è una legenda
plt.legend()
plt.savefig('grafici/temperature.png')

plt.clf()

# Dati per il forward (asse x traslato di -7)
x_forwardK = dataforwardK[:, 0] - 7
y_forwardK = dataforwardK[:, 1]

# Dati per il backward (asse x traslato di -1)
x_backwardK = databackwardK[:, 0] - 1
y_backwardK = databackwardK[:, 1]

# Plot 
plt.plot(x_forwardK, y_forwardK, label="Forward simulation", color="blue", marker='o', markersize=3)
plt.plot(x_backwardK, y_backwardK, label="Backward simulation", color="orange", marker='o', markersize=3)
plt.xlabel("Direzione del tempo")
plt.ylabel("Energia cinetica")
plt.title('Verifica della Reversibilità Temporale:\nAndata e Ritorno dell\'energia cinetica' )
plt.grid(True)
plt.text(0.59, 0.7,"Nota: in 0 finisce la \nsimulazione forward e inizia\nla simulazione backward", transform=plt.gca().transAxes, fontsize=9, fontweight='bold')#aggiungo un commento al grafico, non è una legenda
plt.legend()
plt.savefig('grafici/energia cinetica.png')

plt.clf()

# Dati per il forward (asse x traslato di -7)
x_forwardP = dataforwardP[:, 0] - 7
y_forwardP = dataforwardP[:, 1]

# Dati per il backward (asse x traslato di -1)
x_backwardP = databackwardP[:, 0] - 1
y_backwardP = databackwardP[:, 1]

# Plot 
plt.plot(x_forwardP, y_forwardP, label="Forward simulation", color="blue", marker='o', markersize=3)
plt.plot(x_backwardP, y_backwardP, label="Backward simulation", color="orange", marker='o', markersize=3)
plt.xlabel("Direzione del tempo")
plt.ylabel("Energia potenziale")
plt.title('Verifica della Reversibilità Temporale:\nAndata e Ritorno dell\'energia potenziale' )
plt.grid(True)
plt.text(0.59, 0.2,"Nota: in 0 finisce la \nsimulazione forward e inizia\nla simulazione backward", transform=plt.gca().transAxes, fontsize=9, fontweight='bold')#aggiungo un commento al grafico, non è una legenda
plt.legend()
plt.savefig('grafici/energia potenziale.png')
