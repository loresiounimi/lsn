import matplotlib.pyplot as plt
import numpy as np

#carico i dati dal file
f1, x1, error1 = np.loadtxt("<r>.txt", delimiter=' ', unpack=True) # unpack fa sì che np.loadtxt restituisca le colonne come variabili separate invece che come un unico array 2D.
#creo il grafico
plt.figure()
#errorbar traccia i punti (x,y) con barre di errore yerr. label è l'etichetto per la legenda
plt.errorbar(x1,f1-1/2,yerr=error1)
plt.xlabel('throws')
plt.ylabel('<r>-1/2')
plt.title('<r>-1/2')
plt.grid(True)
plt.savefig("grafici/<r>.png")

f2, x2, error2 = np.loadtxt("r_sigma.txt", delimiter=' ', unpack=True)
plt.figure()
plt.errorbar(x2,f2-1/12,yerr=error2)
plt.xlabel('throws')
plt.ylabel('sigma^2-1/12')
plt.title(r'$\sigma^2$-1/12')
plt.grid(True)
plt.savefig("grafici/r_sigma.png")

