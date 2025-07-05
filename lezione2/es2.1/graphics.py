import matplotlib.pyplot as plt
import numpy as np

#carico i dati dal file
f1, x1, err1 = np.loadtxt("intmedia.txt", delimiter=' ', unpack=True)
#creo il grafico
plt.figure()
#errorbar traccia i punti (x,y) con barre di errore yerr. label è l'etichetto per la legenda
plt.errorbar(x1,f1-1,yerr=err1, label='Campionamento uniforme')
#r serve per scrivere il testo in Latex
plt.title(r'Stima di $I = \int_0^1 \frac{\pi}{2} \cos\left(\frac{\pi x}{2}\right) dx $ con campionamento uniforme')
plt.xlabel('Lanci')
plt.ylabel('I-1')
plt.legend()
plt.grid(True)
plt.savefig('grafici/media.png')

plt.close()

f2, x2, err2 = np.loadtxt("intsampling.txt", delimiter=' ', unpack=True)
plt.figure()
plt.errorbar(x2,f2-1,yerr=err2, label='Importance sampling')
plt.title(r'Stima di $I = \int_0^1 \frac{\pi}{2} \cos\left(\frac{\pi x}{2}\right) dx $ con importance sampling')
plt.xlabel('Lanci')
plt.ylabel('I-1')
plt.legend()
plt.grid(True)
plt.savefig('grafici/sampling.png')