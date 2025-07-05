import matplotlib.pyplot as plt
import numpy as np

#chiquadro singoli
#carico i dati
chi = np.loadtxt("chiquadro_singoli.txt", delimiter=' ', unpack=True)
#creo il grafico
plt.figure()
#creo l'istogramma riempiendolo coi dati
plt.hist(chi, bins=30,  color='skyblue', edgecolor='black')
plt.xlabel(r'$\chi^2$')
plt.ylabel('Occorrenze')
plt.title(r'Chi quadro singoli: $\chi^2$')
plt.grid(True)
plt.savefig("grafici/chiquadro_singoli.png")

#chiquadro con medie cumulative e data blocking
f1, error, x1 = np.loadtxt("chiquadro.txt", delimiter=' ', unpack=True)
plt.figure()
plt.errorbar(x1,f1-100,yerr=error, linestyle='-')
plt.xlabel('Step')
plt.ylabel(r'$\chi^2$-100')
plt.title(r'Chi quadro: $\chi^2$')
plt.grid(True)
plt.savefig("grafici/chiquadro.png")
