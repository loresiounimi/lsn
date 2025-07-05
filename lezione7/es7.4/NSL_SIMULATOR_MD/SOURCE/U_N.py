import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data =np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")
block = data[:,0]
U = data[:,2]
error = data[:,3]
#creo il grafico
plt.figure(figsize=(8, 5))
#plotto i dati con gli errori
plt.errorbar(block, U, yerr=error, label="U/N simulazione MD")

plt.xlabel("Blocco")
plt.ylabel("U/N")
plt.title("U/N simulazione MD")
plt.grid(True)
plt.legend()
plt.savefig('grafici/U_N')

