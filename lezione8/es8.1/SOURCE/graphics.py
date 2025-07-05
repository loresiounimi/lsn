import numpy as np
import matplotlib.pyplot as plt

#calcolo della soluzione matriciale di psi(x)
# Funzione potenziale V(x)
def Vpot(x):
    return (x**2 - 2.5)*x**2  

# Costanti fisiche (unità naturali: ħ = 1, m = 1)
hbar = 1
m = 1
a = 10          # Intervallo totale in cui risolviamo il problema ([-a/2, a/2])
N = 1000        # Numero di punti nella discretizzazione

# Creazione del dominio spaziale x e del passo dx
x = np.linspace(-a/2, a/2, N)
dx = x[1] - x[0]            # Passo spaziale
V = Vpot(x)                 # Valori del potenziale sul dominio x

# Costruzione della matrice delle differenze finite centrali per la derivata seconda:
# f'' ≈ (f_{i+1} - 2*f_i + f_{i-1}) / dx^2
CDiff = (
    np.diag(np.ones(N-1), -1)  # diagonale sotto la principale
    - 2 * np.diag(np.ones(N), 0)  # diagonale principale
    + np.diag(np.ones(N-1), 1)  # diagonale sopra la principale
)

# Costruzione della matrice hamiltoniana:
# H = -ħ²/(2m) * d²/dx² + V(x)
H = (-(hbar**2) * CDiff) / (2 * m * dx**2) + np.diag(V)

# Calcolo degli autovalori (energie) e autovettori (funzioni d'onda) dell'Hamiltoniano
E, psi = np.linalg.eigh(H)

# Trasposizione della matrice psi: ogni riga sarà una funzione d'onda 
# Normalizzazione delle funzioni d'onda (integrale numerico ≈ somma dei quadrati * dx = 1)
psi = np.transpose(psi)
psi = psi / np.sqrt(dx)

# x_samples: array di configurazioni campionate
x_samples = np.loadtxt("../OUTPUT/x_sampled.dat")  

#migliori parametri di mu e sigma
data = np.loadtxt("../../es8.2/OUTPUT/best_parameters.dat", comments="#")
mu = data[0]
sigma = data[1]
#calcolo di psi normalizzata con tali parametri
x1 = np.linspace(-3, 3, 1000)
psi_T = np.exp(-(x1 - mu)**2 / (2 * sigma**2)) + np.exp(-(x1 + mu)**2 / (2 * sigma**2))
prob_density = psi_T**2
prob_density /= np.trapz(prob_density, x1)  # normalizza l'area sotto la curva


# Plotto i dati
plt.figure(figsize=(8,5))
scale = 0.3
plt.plot(x, scale*V, color="Black", label="Potenziale") # potenziale
plt.plot(x,(psi[0])**2, label= "Soluzione con equazione matriciale")#soluazione equazione matriciale
plt.hist(x_samples, bins=100, density=True, alpha=0.6, label='Campionata $|\Psi_T(x)|^2$')# Istogramma normalizzato della distribuzione campionata
plt.plot(x1, prob_density, label='Analitica $|\Psi_T(x)|^2$', color='red')#soluzione analitica con mu e sigma best
plt.title(r"$|\Psi_T(x)|^2$")
plt.xlabel("x")
plt.grid(True)
#imposto i limiti del grafico
plt.xlim((-3,3))
plt.ylim((-0.6,1.1))
plt.tight_layout()
plt.legend()
plt.savefig('grafici/psiquadro.png')