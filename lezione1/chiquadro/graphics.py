import matplotlib.pyplot as plt
import numpy as np

f, x = np.loadtxt("chiquadro_singoli.txt", delimiter=' ', unpack=True)
plt.figure()
plt.plot(x,f-100,linestyle='-')
plt.xlabel('Step')
plt.ylabel(r'$\chi^2$-100')
plt.title(r'Chi quadro singoli: $\chi^2$')
plt.grid(True)
plt.savefig("grafici/chiquadro_singoli.png")

f1, error, x1 = np.loadtxt("chiquadro.txt", delimiter=' ', unpack=True)
plt.figure()
plt.errorbar(x1,f1-100,yerr=error, linestyle='-')
plt.xlabel('Step')
plt.ylabel(r'$\chi^2$-100')
plt.title(r'Chi quadro: $\chi^2$')
plt.grid(True)
plt.savefig("grafici/chiquadro.png")
