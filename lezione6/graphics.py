import numpy as np
import matplotlib.pyplot as plt
#parametri curve teoriche
# Numero di punti in cui valutare la curva teorica 
points = 100
# Temperature tra 0.2 e 3.0 
T = np.linspace(0.2, 3.0, num=points)
# Inversa della temperatura: beta = 1/T
beta = 1 / T
# Costante di accoppiamento tra spin (J > 0 = ferromagnetico)
J = 1.0
# Numero di spin nel sistema (serve per le correzioni finite-size)
Ns = 50
# Calcolo di tanh(J / T), funzione iperbolica ricorrente nella soluzione esatta del modello di Ising 1D
th = np.tanh(J / T)
# th elevato alla Ns-esima potenza: tanh(J/T)^N 
thN = th ** Ns
# Calcolo di coth(J/T) = 1 / tanh(J/T)
ch = 1 / th

#carico tutti i dati da entrambe le simulazioni
data_m = np.loadtxt("NSL_SIMULATOR_metropolis/OUTPUT/total_energy.dat", comments="#")
data1_m = np.loadtxt("NSL_SIMULATOR_metropolis/OUTPUT/specific_heat.dat", comments="#")
data2_m = np.loadtxt("NSL_SIMULATOR_metropolis/OUTPUT/susceptibility.dat", comments="#")
data3_m = np.loadtxt("NSL_SIMULATOR_metropolis/OUTPUT/magnetization.dat", comments="#")
data_g = np.loadtxt("NSL_SIMULATOR_gibbs/OUTPUT/total_energy.dat", comments="#")
data1_g = np.loadtxt("NSL_SIMULATOR_gibbs/OUTPUT/specific_heat.dat", comments="#")
data2_g = np.loadtxt("NSL_SIMULATOR_gibbs/OUTPUT/susceptibility.dat", comments="#")
data3_g = np.loadtxt("NSL_SIMULATOR_gibbs/OUTPUT/magnetization.dat", comments="#")
#preparo gli array che li conterranno
x = []
y_internal_energy_m = []
y_internal_energy_g = []
error_internal_energy_m = []
error_internal_energy_g = []
y_specific_heat_m = []
y_specific_heat_g = []
error_specific_heat_m = []
error_specific_heat_g = []
y_chi_m = []
y_chi_g = []
error_chi_m = []
error_chi_g = []
y_M_m =[]
y_M_g =[]
error_M_m =[]
error_M_g =[]
#ciclo sulle temperature per riempire gli array. prendo l'ultima media cumulativa e l'ultima incertezza dei 20 blocchi di ogni simulazione a diversa T
for i in range(26):
    x.append(0.5+i*(0.1))#temperature della simulazione
    y_internal_energy_m.append(data_m[19+20*(i),2])
    y_internal_energy_g.append(data_g[19+20*(i),2])
    error_internal_energy_m.append(data_m[19+20*(i),3])
    error_internal_energy_g.append(data_g[19+20*(i),3])
    y_specific_heat_m.append(data1_m[19+20*(i),2])
    y_specific_heat_g.append(data1_g[19+20*(i),2])
    error_specific_heat_m.append(data1_m[19+20*(i),3])
    error_specific_heat_g.append(data1_g[19+20*(i),3])
    y_chi_m.append(data2_m[19+20*(i),2])
    y_chi_g.append(data2_g[19+20*(i),2])
    error_chi_m.append(data2_m[19+20*(i),3])
    error_chi_g.append(data2_g[19+20*(i),3])
    y_M_m.append(data3_m[19+20*(i),2])
    y_M_g.append(data3_g[19+20*(i),2])
    error_M_m.append(data3_m[19+20*(i),3])
    error_M_g.append(data3_g[19+20*(i),3])


# Formula dell'energia media per spin del modello di Ising 1D con N spin 
e = -J * (th + ch * thN) / (1 + thN)

#plotto l'energia interna per particella teorica e quella ottenuta con le simulazioni
plt.errorbar(x, y_internal_energy_m, yerr=error_internal_energy_m, fmt='-', label="Simulated metropolis U/N", color='red', ecolor='black', capsize=3)
plt.errorbar(x, y_internal_energy_g, yerr=error_internal_energy_g, fmt='-', label="Simulated gibbs U/N", color='green', ecolor='yellow', capsize=3)
plt.plot(T, e, label="Theoretical U/N", color='blue')
plt.title('Ising 1D, internal energy')
plt.xlabel('T')
plt.ylabel('U/N')
plt.legend()
plt.savefig('grafici/internal_energy')

plt.clf()

#cv in funzione di T teorico
heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)
#plotto cv ottenuto con simulazioni e cv teorico
plt.errorbar(x, y_specific_heat_m, yerr=error_specific_heat_m, fmt='-', label="Simulated metropolis cv", color='red', ecolor='black', capsize=3)
plt.errorbar(x, y_specific_heat_g, yerr=error_specific_heat_g, fmt='-', label="Simulated gibbs cv", color='green', ecolor='yellow', capsize=3)
plt.plot(T, heat, label="Theoretical cv", color='blue')
plt.title('Ising 1D, Heat capacity')
plt.xlabel('T')
plt.ylabel('C')
plt.legend()
plt.savefig('grafici/heat_capacity')
# Nota che cv è la varianza dell'energia. A basse T le transizioni energetiche sono rare, quindi le configurazioni esplorate sono molto simili e poche e quindi le fluttuazioni di E sono poco rappresentative e sottostimate -> a basse T cv simulato si discosta da cv teorico e risulta sottostimato, essendo sottostimate le fluttuazioni di E e quindi anche la sua varianza che è legata a cv.
plt.clf()

#suscettività in funzione di T teorica
X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
#plotto suscettività ottenuta con simulazioni e teorica
plt.errorbar(x, y_chi_m, yerr=error_chi_m, fmt='-', label="Simulated metropolis chi", color='red', ecolor='black', capsize=3)
plt.errorbar(x, y_chi_g, yerr=error_chi_g, fmt='-', label="Simulated gibbs chi", color='green', ecolor='yellow', capsize=3)
plt.plot(T, X, label="Theoretical chi", color='blue')
plt.title('Ising 1D, Susceptibility')
plt.xlabel('T')
plt.ylabel('$\chi$')
plt.legend()
plt.savefig('grafici/susceptibility')

plt.clf()

#magnetizzazione in funzione di T teorica
h=0.02 #campo esterno

l1 = np.exp(beta*J)*np.cosh(beta*h)+np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
l2 = np.exp(beta*J)*np.cosh(beta*h)-np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
Z = l1**Ns + l2**Ns
M = (np.exp(beta*J)*np.sinh(beta*h)*((l1**(Ns-1))*(1+np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))) 
        + (l2**(Ns-1))*(1-np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J)))))/(Z)

#plotto magnetizzazione ottenuta con simulazioni e teorica
plt.errorbar(x, y_M_m, yerr=error_M_m, fmt='-', label="Simulated metropolis M", color='red', ecolor='black', capsize=3)
plt.errorbar(x, y_M_g, yerr=error_M_g, fmt='-', label="Simulated gibbs M", color='green', ecolor='yellow', capsize=3)
plt.plot(T, M, label="Theoretical M", color='blue')
plt.title('Ising 1D, magnetization M with h = 0.02')
plt.xlabel('T')
plt.ylabel('$M$')
plt.legend()
plt.savefig('grafici/magnetization')