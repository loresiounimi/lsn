import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussiana(x, A, mu, sigma):
  return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
def lorentziana(x, A, mu, gamma):
    return A * (1/np.pi) * (gamma / ((x - mu)**2 + gamma**2))

a = np.array([1, 1, 8, 97])
std = np.loadtxt("standardisto.txt")
for i in range(4):
  plt.figure()
  if(i==3):
   counts, bins = np.histogram(std[:,i], bins=100, density='true')
   bin_centers = (bins[:-1] + bins[1:]) / 2
   p0 = [np.max(counts), np.mean(std[:,i]), np.std(std[:,i])] #guess iniziale dei parametri della gaussiana A ~ max(counts), mu ~ meanstd[:,3]), sigma ~ std([:,3])
   popt, _ = curve_fit(gaussiana, bin_centers, counts, p0=p0) #popt è array con i parametri ottimali trovati dal fit. _ indica che ignoro pcov e cioè la matrice di covarianza dei parametri del fit, serve per l'incertezza
  plt.hist(std[:, i], bins=100, color='skyblue', edgecolor='black', density='true')  # Istogramma
  if(i==3):
    x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000) #creo 1000 punti tra il centro del primo (0) e ultimo (-1) bin
    plt.plot(x_fit, gaussiana(x_fit, *popt), 'r-', label='Fit Gaussiana') #plotto la gaussiana fittata
    plt.legend()
  plt.ylabel('Conteggi')  # Etichetta asse x
  plt.xlabel(f'Sn N={i+a[i]}')  # Etichetta asse y
  plt.title(f'Istogramma Sn N={i+a[i]}')  # Titolo del grafico   
  plt.savefig(f'grafici/standard_isto N={i+a[i]}')
  plt.close()

exp = np.loadtxt("expisto.txt")
for i in range(4):
 plt.figure() 
 if(i==3):
   counts, bins = np.histogram(exp[:,i], bins=100, density='true')
   bin_centers = (bins[:-1] + bins[1:]) / 2
   p0 = [np.max(counts), np.mean(exp[:,i]), np.std(exp[:,i])] 
   popt, _ = curve_fit(gaussiana, bin_centers, counts, p0=p0) 
 plt.hist(exp[:, i], bins=100, color='skyblue', edgecolor='black', density='true')  
 if(i==3):
    x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000) 
    plt.plot(x_fit, gaussiana(x_fit, *popt), 'r-', label='Fit Gaussiana') 
    plt.legend()
 plt.ylabel('Conteggi')  
 plt.xlabel(f'Sn N={i+a[i]}')  
 plt.title(f'Istogramma Sn N={i+a[i]}')  
 plt.savefig(f'grafici/exp_isto N={i+a[i]}')
 plt.close()

cauchy = np.loadtxt("cauchyisto.txt")
for i in range(4):
 plt.figure()  
 if(i==3):
   counts, bins = np.histogram(cauchy[:,i], bins=100, range=(-25,25), density='true')
   bin_centers = (bins[:-1] + bins[1:]) / 2
   p0 = [np.max(counts), np.median(cauchy[:, i]), (np.percentile(cauchy[:, i], 75) - np.percentile(cauchy[:, i], 25)) / 2] #np.percentile calcola l'ennesimo percentile, cioè il valore della x sotto il quale sta il N% dei dati. la stima di gamma fatta è l'IQR ( cioè l’intervallo che contiene il 50% centrale dei dati) diviso 2 
   popt, _ = curve_fit(lorentziana, bin_centers, counts, p0=p0)
 plt.hist(cauchy[:, i], range=(-25,25), bins=100, color='skyblue', edgecolor='black', density='true')  
 if(i==3):
    x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000) 
    plt.plot(x_fit, lorentziana(x_fit, *popt), 'r-', label='Fit Lorentziana') 
    plt.legend()
 plt.ylabel('Conteggi')  
 plt.xlabel(f'Sn N={i+a[i]}') 
 plt.title(f'Istogramma Sn N={i+a[i]}')  
 plt.savefig(f'grafici/cauchy_isto N={i+a[i]}')
 plt.close()

