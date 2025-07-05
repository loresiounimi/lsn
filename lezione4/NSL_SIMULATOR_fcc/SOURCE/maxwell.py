import numpy as np
import matplotlib.pyplot as plt

# parametri noti
T_star = 1.34651 # temperatura usata per il fit. è stata ottenuta prendendo dati dopo un iniziale loop di termalizzazione
num_blocks = 20   # numero di blocchi in cui è suddivisa la simulazione
bins_per_block = 30  # numero di bin (intervalli di velocità) per ciascun blocco

# carico i dati dal file pofv.dat 
data = np.loadtxt("../OUTPUT/pofv.dat", comments="#")

# larghezza del bin: il primo valore della prima riga indica il centro del primo bin,
# quindi bin_width è il doppio (per coprire tutta la larghezza del bin). è uguale per ogni bin
bin_width = data[0, 0] * 2

# definizione della distribuzione teorica di Maxwell-Boltzmann
def maxwell_boltzmann(v, T_star):
    prefactor = 4 * np.pi * v**2 / ((2 * np.pi * T_star) ** (3/2))
    return prefactor * np.exp(-v**2 / (2 * T_star))

# ciclo per generare e salvare un grafico per ciascun blocco
for i in range(num_blocks):
    # estraggo i dati del blocco corrente (30 righe alla volta)
    block_data = data[i * bins_per_block : (i + 1) * bins_per_block]

    # prendo i centri dei bin, le altezze degli istogrammi e gli errori
    bin_centers = block_data[:, 0]  # valori centrali delle velocità
    heights = block_data[:, 2]      # frequenze (probabilità) osservate
    errors = block_data[:, 3]       # incertezze sulle frequenze

    # velocità da usare per il calcolo della curva teorica
    v = bin_centers
    theoretical = maxwell_boltzmann(v, T_star)

    # creo una nuova figura
    plt.figure()
    
    # disegno l'istogramma con le barre (con errori)
    plt.bar(bin_centers, heights, width=bin_width, yerr=errors, capsize=3,
            color='skyblue', edgecolor='black', label='Simulazione')

    # disegno la curva teorica di Maxwell-Boltzmann
    plt.plot(v, theoretical, 'r-', linewidth=2, label='Maxwell-Boltzmann teorica')

    # etichette degli assi e titolo
    plt.xlabel('Velocità')
    plt.ylabel('Probabilità')
    plt.title(f'Distribuzione velocità progressiva - Blocco {i+1}')
    
    # griglia e legenda
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # salva la figura con nome progressivo (es. maxwell_block_01.png, ..._02.png, ...)
    filename = f"grafici/maxwell_block_{i+1:02}.png"  
    plt.savefig(filename)
    plt.close()  

#carico i dati dell'energia cinetica
data = np.loadtxt("../OUTPUT/kinetic_energy.dat", comments="#")
blocco = data[:,0] 
K_blocco = data[:,1] 
K_cumulativa = data[:,2] 
K_errore = data[:,3] 

# creo il grafico
plt.figure()

# energia cinetica media per blocco
plt.plot(blocco, K_blocco, label="K blocco", color='skyblue')

# media cumulativa
plt.errorbar(blocco, K_cumulativa, yerr=K_errore, color='red',fmt='o-',capsize=4, label="K cumulativa ± errore")

# etichette e titolo
plt.xlabel("Blocco")
plt.ylabel("Energia cinetica per particella")
plt.title("Energia cinetica: media per blocco e cumulativa")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('grafici/energia_cinetica.png')

plt.close()

#carico i dati dell'energia potenziale
data = np.loadtxt("../OUTPUT/potential_energy.dat", comments="#")
blocco = data[:,0] 
U_blocco = data[:,1] 
U_cumulativa = data[:,2] 
U_errore = data[:,3] 

# creo il grafico
plt.figure()

# energia potenziale media per blocco
plt.plot(blocco, U_blocco, label="K blocco", color='skyblue')

# media cumulativa
plt.errorbar(blocco, U_cumulativa, yerr=U_errore, color='red',fmt='o-',capsize=4, label="U cumulativa ± errore")

# etichette e titolo
plt.xlabel("Blocco")
plt.ylabel("Energia potenziale per particella")
plt.title("Energia potenziale: media per blocco e cumulativa")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('grafici/energia_potenziale.png')

plt.close()

#carico i dati dell'energia totale
data = np.loadtxt("../OUTPUT/total_energy.dat", comments="#")
blocco = data[:,0] 
T_blocco = data[:,1] 
T_cumulativa = data[:,2] 
T_errore = data[:,3] 

# creo il grafico
plt.figure()

# energia totale media per blocco
plt.plot(blocco, T_blocco, label="T blocco", color='skyblue')

# media cumulativa
plt.errorbar(blocco, T_cumulativa, yerr=T_errore, color='red',fmt='o-',capsize=4, label="T cumulativa ± errore")

# etichette e titolo
plt.xlabel("Blocco")
plt.ylabel("Energia totale per particella")
plt.title("Energia totale: media per blocco e cumulativa")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('grafici/energia_totale.png')

plt.close()

#carico i dati della temperatura
data = np.loadtxt("../OUTPUT/temperature.dat", comments="#")
blocco = data[:,0] 
temp_blocco = data[:,1] 
temp_cumulativa = data[:,2] 
temp_errore = data[:,3] 

# creo il grafico
plt.figure()

# temperatura media per blocco
plt.plot(blocco, temp_blocco, label="temperature del blocco", color='skyblue')

# media cumulativa
plt.errorbar(blocco, temp_cumulativa, yerr=temp_errore, color='red',fmt='o-',capsize=4, label="temperatura cumulativa ± errore")

# etichette e titolo
plt.xlabel("Blocco")
plt.ylabel("temperatura per particella")
plt.title("Temperatura: media per blocco e cumulativa")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('grafici/temperatura.png')

plt.close()