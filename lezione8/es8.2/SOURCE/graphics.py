import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data =np.loadtxt("../OUTPUT/SAevolution.dat", comments="#")
data1 =np.loadtxt("../OUTPUT/SAparameters.dat", comments="#")
SAstep = data[:,0].astype(int)
H_medio = data[:,1]
H_error = data[:,2]
mu = data1[:,0]
sigma = data1[:,1]

#creo il grafico con i valori di <H> durante la simulazione SA
plt.figure(figsize=(8, 5))
plt.errorbar(SAstep, H_medio, yerr=H_error, label="<H>", fmt='-', capsize=3)
plt.xlabel("SAstep")
plt.ylabel("<H>")
plt.title("<H> in funzione dei passi di SA")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('grafici/SAevolution.png')

plt.clf()

#creo il grafico con l'evoluzione dei parametri mu e sigma
plt.figure(figsize=(8, 6))
plt.plot(mu, sigma, linestyle='-', marker = 'o', color='blue' )
# Evidenzio il punto iniziale (primo)
plt.scatter(mu[0], sigma[0], color='green', s=100, label='Inizio', zorder=5)
# Evidenzio il punto finale (ultimo)
plt.scatter(mu[-1], sigma[-1], color='red', s=100, label='Fine', zorder=5)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\sigma$')
plt.title('Evoluzione dei parametri mu e sigma')
plt.legend()
plt.tight_layout()
plt.savefig('grafici/SAparameters.png')