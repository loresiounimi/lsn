import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data =np.loadtxt("../OUTPUT/pressure.dat", comments="#")
block = data[:,0]
P = data[:,2]
error = data[:,3]
#creo il grafico
plt.figure(figsize=(8, 5))
#plotto i dati con gli errori
plt.errorbar(block, P, yerr=error, label="P simulazione MD")

plt.xlabel("Blocco")
plt.ylabel("P")
plt.title("P simulazione MD")
plt.grid(True)
plt.legend()
plt.savefig('grafici/P')