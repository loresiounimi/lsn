import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data = np.loadtxt("../OUTPUT/<H>.dat", comments="#")
Nblocks = data[:,0]
H_medio = data[:,2]
H_error = data[:,3]

#creo il grafico
plt.figure(figsize=(8, 5))
plt.errorbar(Nblocks, H_medio, yerr=H_error, label="<H>", fmt='-o', capsize=3)
plt.xlabel("MCblocks")
plt.ylabel("<H>")
plt.title("<H> with best parameters. MC simulation ")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('grafici/<H>_best_parameters.png')