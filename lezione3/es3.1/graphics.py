import matplotlib.pyplot as plt
import numpy as np

#carico i dati dal file txt
x1, f1_0, err1_0, f1_1, err1_1 = np.loadtxt("callprice.txt", delimiter=' ', unpack=True)
#creo il grafico
plt.figure()
#errorbar traccia i punti (x,y) con barre di errore yerr. label è l'etichetto per la legenda
plt.errorbar(x1,f1_0,yerr=err1_0, label='Wiener 1 passo', color='blue')
plt.errorbar(x1,f1_1,yerr=err1_1, label='Wiener 100 passi', color='red')
#traccio una linea orizzontale sul valore della soluazione analitica di Black-Scholes
plt.axhline(y=14.975790778311286, color='green', linestyle='--', label='Soluzione analitica di Black-Scholes')
plt.xlabel('lanci')
plt.ylabel('<C[S(0),0]>')
plt.title('Prezzo dell\'opzione call in funzione del numero di lanci')
plt.legend()
plt.grid(True)

plt.savefig('grafici/callprice.png')

x2, f2_0, err2_0, f2_1, err2_1 = np.loadtxt("putprice.txt", delimiter=' ', unpack=True)

plt.figure()
plt.errorbar(x2,f2_0,yerr=err2_0, label='Wiener 1 passo', color='blue')
plt.errorbar(x2,f2_1,yerr=err2_1, label='Wiener 100 passi', color='red')
plt.axhline(y=5.4595325819072364, color='green', linestyle='--', label='Soluzione analitica di Black-Scholes')
plt.xlabel('lanci')
plt.ylabel('<P[S(0),0]>')
plt.title('Prezzo dell\'opzione put in funzione del numero di lanci')
plt.legend()
plt.grid(True)

plt.savefig('grafici/putprice.png')