import matplotlib.pyplot as plt
import numpy as np
import math

#carico i dati
x1, f1, error1 = np.loadtxt("Buffon.txt", delimiter=' ', unpack=True) 
#creo il grafico
plt.figure()
plt.errorbar(x1,f1-math.pi,yerr=error1)
plt.xlabel('throws')
plt.ylabel(r'value-$\pi$')
plt.title(r'$\pi$ from Buffon experiment')
plt.grid(True)
plt.savefig("grafici/Buffon.png")