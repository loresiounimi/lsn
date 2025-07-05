import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#definisco la funzione da fittare f(x)=k*(radice di N)+c
def funzione_fit(x, k, c):
    return (k * np.sqrt(x))+c

#carico i dati dal file
x1, f1, err1 = np.loadtxt("randomdiscretewalk.txt", delimiter=' ', unpack=True)
#faccio il fit della funzione sui dati del file
popt,_ = curve_fit(funzione_fit, x1, f1)# l'underscore serve a python per ignorare la covarianza, cioè il calcolo dell'errore sul fit, che a me in questo caso non interessa
k, c = popt #setto i parametri della funzione con quelli del fit effettuato

x_fit = np.linspace(min(x1), max(x1), 1000) #numero e valore dei punti su cui eseguire il fit
y_fit = funzione_fit(x_fit, k, c) #fit della funzione
#creo il rafico
plt.figure()
plt.errorbar(x1,f1,yerr=err1, label='Dati')
#errorbar traccia i punti (x,y) con barre di errore yerr. label è l'etichetto per la legenda
plt.plot(x_fit, y_fit, label=rf'Fit: $k\sqrt{{N}} + c$, $k = {k:.3e}, $c = {c:.3f} ', color='red') #la notazione {k:.3e} mostra k in notazione scientifica con 3 cifre decimali, la notazione {c:.3f} mostra c con 3 cifre decimali. $ mostra la scritta in corsivo
#r serve per scrivere il testo in Latex
plt.xlabel('Numero di passi N')
plt.ylabel(r'$\sqrt{\langle |\vec{r}_N|^2 \rangle}$')
plt.title('Random walk 3D su reticolo cubico')
plt.legend()
plt.grid(True)

plt.savefig('grafici/randomdiscretewalk.png')

plt.close()

x2, f2, err2 = np.loadtxt("randomcontinuumwalk.txt", delimiter=' ', unpack=True)

popt, _ = curve_fit(funzione_fit, x2, f2)
k, c = popt

x_fit = np.linspace(min(x2), max(x2), 1000) 
y_fit = funzione_fit(x_fit, k, c) 
plt.figure()
plt.errorbar(x2,f2,yerr=err2, label='Dati')
plt.plot(x_fit, y_fit, label=f'Fit: $k\sqrt{{N}} + c$, $k = {k:.3e}, $c = {c:.3f} ', color='red')
plt.xlabel('Numero di passi N')
plt.ylabel(r'$\sqrt{\langle |\vec{r}_N|^2 \rangle}$')
plt.title('Random walk 3D continuo')
plt.legend()
plt.grid(True)

plt.savefig('grafici/randomcontinuumwalk.png')
