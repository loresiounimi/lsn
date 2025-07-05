import numpy as np
import matplotlib.pyplot as plt

#carico i dati
data =np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")
#salvo solo i dati della prima colonna del file, cioè le misure sui singoli step. infatti calcolo l'incertezza in questo codice
un_istantaneo = data[:, 1]
#numero totale degli step
M = len(un_istantaneo)  
#preparo gli array per errore e lunghezza del blocco
error = []
range_blk = []
#ciclo su L lunghezza blocco
for L in range(10, 5001, 100):
    #salvo il valore di L
    range_blk.append(L)
    #numero di steps per blocco
    N = M // L  
    #array per la media di blocco
    mean = []
    #calcolo medie di blocco sui N blocchi
    for i in range (N): 
        mean.append(np.mean(un_istantaneo[L*i:L*(i+1)]))
    #calcolo la media cumulativa delle medie di blocco
    mean = np.array(mean)
    #calcolo dell'errore con data blocking
    error.append(np.sqrt((np.mean(mean**2)-np.mean(mean)**2)/N))
#salvo lunghezza dei blocchi ed incertezze negli array
error = np.array(error)
range_blk = np.array(range_blk)

#creo il grafico
plt.figure(figsize=(8, 5))
#plotto gli errori in funzione della lunghezza dei blocchi in scala logaritmica
plt.plot(range_blk, error, marker='o', label="Statistical uncertainties")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Larghezza del blocco")
plt.ylabel("Errore statistico")
plt.title("Errore statistico in funzione della larghezza del blocco")
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.legend()
plt.savefig('grafici/errorestatistico')
#all'inizio con L l'errore è sottostimato a causa del gran numero di blocchi che dunque sottostima l'errore nel calcolo. Inoltre i blocchi sono troppo corti per rompere le correlazioni temporali tra i dati e pertanto le medie di blocco non sono indipendenti; anche questo contribuisce alla sottostima. poi l'errore si stabilizza attorno a un certo valore oscillando. sono anch'esse fluttuazioni dovute alla natura casuale dei valori di U/N. 
#se aumentassi ancora di più L fino a diventare paragonabile a M troverei fluttuazioni sempre più "violente" perchè il numero di blocchi diventa molto piccolo e l'errore di conseguenza sovrastimato e molto sensibile alla natura casuale dei valori di U/N. In pratica l'incertezza sull'errore diventa significativamente grande.