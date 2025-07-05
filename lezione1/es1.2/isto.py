import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Funzione gaussiana per il fit
def gaussiana(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Funzione lorentziana (Cauchy) per il fit
def lorentziana(x, A, mu, gamma):
    return A * (1/np.pi) * (gamma / ((x - mu)**2 + gamma**2))

# Array con valori per etichettare i grafici: N = i + a[i]
a = np.array([1, 1, 8, 97])

# Caricamento dei dati 
std = np.loadtxt("standardisto.txt")

# Loop sui 4 histogrammi da plottare
for i in range(4):
    plt.figure()  

    # Solo per l'ultima colonna (i == 3) faccio anche il fit
    if i == 3:
        # Calcolo l'istogramma normalizzato (density=True)da cui prendere i dati per il fit
        counts, bins = np.histogram(std[:, i], bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2  # Centri dei bin per il fit presi a metà di ogni bin

        # Parametri iniziali per il fit:
        # - ampiezza (A): stimata come il valore massimo dell'istogramma (np.max(counts))
        # - media (mu): stimata come la media dei dati nella colonna (np.mean(std[:, i]))
        # - deviazione standard (sigma): stimata come la deviazione standard dei dati (np.std(std[:, i]))
        p0 = [np.max(counts), np.mean(std[:, i]), np.std(std[:, i])]

        # Fit dei dati con curva gaussiana
        #popt è array con i parametri ottimali trovati dal fit. _ indica che ignoro pcov e cioè la matrice di covarianza dei parametri del fit, serve per l'incertezza
        popt, _ = curve_fit(gaussiana, bin_centers, counts, p0=p0)

    # Istogramma dei dati normalizzato (density=true)
    plt.hist(std[:, i], bins=100, color='skyblue', edgecolor='black', density=True)

    # Plot del fit (solo per i == 3)
    if i == 3:
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000) #creo 1000 punti tra il centro del primo (0) e ultimo (-1) bin
        plt.plot(x_fit, gaussiana(x_fit, *popt), 'r-', label='Fit Gaussiana')
        plt.legend()

    # Etichette e salvataggio
    plt.ylabel('Conteggi')
    plt.xlabel(f'Sn N={i + a[i]}')
    plt.title(f'Istogramma Sn N={i + a[i]}')
    plt.savefig(f'grafici/standard_isto N={i + a[i]}')
    plt.close()  # Chiudo la figura per liberare memoria

# Caricamento dei dati 
exp = np.loadtxt("expisto.txt")

# Stesso procedimento per i dati exp
for i in range(4):
    plt.figure()
    
    if i == 3:
        counts, bins = np.histogram(exp[:, i], bins=100, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        p0 = [np.max(counts), np.mean(exp[:, i]), np.std(exp[:, i])]
        popt, _ = curve_fit(gaussiana, bin_centers, counts, p0=p0)

    plt.hist(exp[:, i], bins=100, color='skyblue', edgecolor='black', density=True)

    if i == 3:
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        plt.plot(x_fit, gaussiana(x_fit, *popt), 'r-', label='Fit Gaussiana')
        plt.legend()

    plt.ylabel('Conteggi')
    plt.xlabel(f'Sn N={i + a[i]}')
    plt.title(f'Istogramma Sn N={i + a[i]}')
    plt.savefig(f'grafici/exp_isto N={i + a[i]}')
    plt.close()

# Caricamento dei dati
cauchy = np.loadtxt("cauchyisto.txt")

# Stesso ciclo per i dati Cauchy
for i in range(4):
    plt.figure()

    if i == 3:
        # Istogramma con range limitato per evitare code infinite di Cauchy
        counts, bins = np.histogram(cauchy[:, i], bins=100, range=(-25, 25), density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # Parametri iniziali: A, mediana, gamma stimato da IQR/2
        p0 = [np.max(counts),
              np.median(cauchy[:, i]),
              (np.percentile(cauchy[:, i], 75) - np.percentile(cauchy[:, i], 25)) / 2]

        # Fit con la funzione lorentziana
        popt, _ = curve_fit(lorentziana, bin_centers, counts, p0=p0)

    # Istogramma con range limitato (evita code troppo estese della Cauchy)
    plt.hist(cauchy[:, i], bins=100, range=(-25, 25), color='skyblue', edgecolor='black', density=True)

    if i == 3:
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        plt.plot(x_fit, lorentziana(x_fit, *popt), 'r-', label='Fit Lorentziana')
        plt.legend()

    plt.ylabel('Conteggi')
    plt.xlabel(f'Sn N={i + a[i]}')
    plt.title(f'Istogramma Sn N={i + a[i]}')
    plt.savefig(f'grafici/cauchy_isto N={i + a[i]}')
    plt.close()
