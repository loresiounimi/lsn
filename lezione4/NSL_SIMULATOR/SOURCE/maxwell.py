import numpy as np
import matplotlib.pyplot as plt

# Carico il file
data = np.loadtxt("../OUTPUT/maxwell.dat", comments="#") #maxwell.dat contiene i dati relativi all'ultimo blocco della simulazione

# Estraggo le colonne
bin_centers = data[:, 0]
heights = data[:, 2]
errors = data[:, 3]

bin_width = data[0,0]*2  

# Creo l’istogramma con barre d’errore
plt.bar(bin_centers, heights, width=bin_width, align='center', yerr=errors, capsize=3, color='skyblue', edgecolor='black', label='simulazione')#align serve per allineamento centrale

# Parametro della temperatura in unità ridotte T
T_star = 0.561875 #preso da temperature.dat, ultima riga

# Calcolo la distribuzione teorica di Maxwell-Boltzmann normalizzata prendendo come ascisse i centri dei bins
v = bin_centers
prefactor = 4 * np.pi * v**2 / ((2 * np.pi * T_star) ** (3/2))
exp_factor = np.exp(-v**2 / (2 * T_star))
theoretical = prefactor * exp_factor

# Plot della funzione teorica
plt.plot(v, theoretical, 'r-', linewidth=2, label=rf'$p(v^*, T^*) = \frac{{4\pi v^{{*2}}}}{{(2\pi \cdot {T_star:.3f})^{{3/2}}}} e^{{-\frac{{v^{2}}}{{2 \cdot {T_star:.3f}}}}}$')

plt.xlabel('Velocità')
plt.ylabel('Probabilità')
plt.title('Istogramma con barre di errore della distribuzione di velocità')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('grafici/maxwell.png')