import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#carico i dati
data =np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")
#variabile dipendente t della funzione di autocorrelazione
time = data[:, 0].astype(int)
#m(t) della funzione di autocorrelazione
un_istantaneo = data[:, 1]

# Calcolo della funzione di autocorrelazione normalizzata
def chi_autocorrelation(m):
    m = np.array(m)
    tmax = int(np.max(time))

    mean_total = np.mean(m)
    mean_sq_total = np.mean(m**2)
    denom = mean_sq_total - mean_total**2

    chi = []
    for t in range(tmax): #range fino a tmax-1
        m1 = m[:tmax - t] #va da m(0) a m(tmax-t-1)
        m2 = m[t:] #va da m(t) a m(tmax-1) 
        #m1 e m2 hanno entrambi lunghezza tmax-t
        prod_mean = np.mean(m1 * m2)
        mean1 = np.mean(m1)
        mean2 = np.mean(m2)
        num = prod_mean - mean1 * mean2
        chi.append(num / denom)
    return np.array(chi)

chi_t = chi_autocorrelation(un_istantaneo)

#definizione della funzione da fittare con l'autocorrelazione
def exp_fit(t, A, t_c):
    return A * np.exp(-t / t_c)

#fit fino al 200000esimo valore, la coda della funzione infatti subisce forti fluttuazioni che "falsano" il fit
params,_ = curve_fit(exp_fit, time, chi_t, p0=[chi_t[0], 1.0])
A_fit, t_c_fit = params

#creo il grafico
plt.figure(figsize=(8, 5))
#plotto la funzione di autocorrelazione e il fit per i primi 2000 steps
plt.plot(time[:2000], exp_fit(time[:2000], *params),label=fr"Fit esponenziale: $A \cdot e^{{-t/t_c}}$, $t_c = {t_c_fit:.2f}$",color="red",linestyle="--") #* davanti a params estrae separatamente i valori nell'array params
plt.plot(time[:2000], chi_t[:2000], label="χ(t)")

plt.xlabel("Tempo discreto t")
plt.ylabel("χ(t)")
plt.title("Funzione di autocorrelazione di U/N")
plt.grid(True)
plt.legend()
plt.savefig('grafici/autocorrelazione')

# per vedere il seguente effetto è necessario plottare i dati su tutto il dataset completo

#il fatto che la coda della funzione non sia zero ma fluttui è causato dalle fluttuazioni statistiche di U/N(t) e anche dalla natura stessa di x(t).
#Quando t si avvicina a tmax, la quantità di dati disponibili per calcolare la funzione di autocorrelazione diminuisce. Di conseguenza, la media e la varianza delle fluttuazioni misurate diventano meno precise, e questo può portare a fluttuazioni più ampie e statisticamente poco affidabili nella funzione di autocorrelazione.
