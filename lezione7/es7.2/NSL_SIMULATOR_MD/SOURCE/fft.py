import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#carico i dati
data =np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")
un_istantaneo = data[:, 1]
time = data[:, 0].astype(int)

#funzione di autocorrelazione fft
def autocorrelation_fft(x, max_lag=None):
    x = x - np.mean(x)  # centro la serie
    n = len(x)
    if max_lag is None or max_lag > n:
        max_lag = n

    f = np.fft.fft(x, n=2*n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    acf /= acf[0]  # normalizzazione
    return acf[:max_lag]

fft_t = autocorrelation_fft(un_istantaneo, max_lag=len(time))

#definizione della funzione da fittare con l'autocorrelazione
def exp_fit(t, A, t_c):
    return A * np.exp(-t / t_c)

params,_ = curve_fit(exp_fit, time, fft_t, p0=[fft_t[0], 1.0])
A_fit, t_c_fit = params

#creo grafico
plt.figure(figsize=(8, 5))
#plotto l'autocorrelazione e la funzione fittata per i primi 2000 steps
plt.plot(time[:2000],exp_fit(time[:2000], *params),label=fr"Fit esponenziale: $A \cdot e^{{-t/t_c}}$, $t_c = {t_c_fit:.2f}$",color="red",linestyle="--") #* davanti a params estrae separatamente i valori nell'array params
plt.plot(time[:2000], fft_t[:2000], label="χ(t) - FFT", linestyle=":", alpha=0.7)

plt.xlabel("Passi")
plt.ylabel("χ(t)") 
plt.title("Funzione di autocorrelazione di U/N") 
plt.grid(True)
plt.legend()
plt.savefig('grafici/autocorrelazionefft')

