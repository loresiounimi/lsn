import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data =np.loadtxt("../OUTPUT/gofr/gofr_block_20.dat", comments="#")
x = data[:,0]
y = data[:,2]
error = data[:,3]
#creo il grafico
plt.figure(figsize=(8, 5))
#plotto i dati con gli errori
plt.errorbar(x, y, yerr=error, label="g(r) simulazione MC")

plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("g(r) simulazione MC")
plt.grid(True)
plt.legend()
plt.savefig('grafici/g(r)')