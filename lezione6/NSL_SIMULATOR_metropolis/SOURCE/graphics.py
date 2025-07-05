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
# Formula dell'energia media per spin del modello di Ising 1D con N spin 
e = -J * (th + ch * thN) / (1 + thN)

#carico i dati ottenuti con la simulzione
data = np.loadtxt("../OUTPUT/total_energy.dat", comments="#")
data1 = np.loadtxt("../OUTPUT/specific_heat.dat", comments="#")
data2 = np.loadtxt("../OUTPUT/susceptibility.dat", comments="#")
data3 = np.loadtxt("../OUTPUT/magnetization.dat", comments="#")
#preparo gli array che li conterranno
x = []
y_internal_energy = []
error_internal_energy = []
y_specific_heat = []
error_specific_heat = []
y_chi = []
error_chi = []
y_M =[]
error_M =[]
#ciclo sulle temperature per riempire gli array. prendo l'ultima media cumulativa e l'ultima incertezza dei 20 blocchi di ogni simulazione a diversa T
for i in range(26):
    x.append(0.5+i*(0.1))#temperature della simulazione
    y_internal_energy.append(data[19+20*(i),2])
    error_internal_energy.append(data[19+20*(i),3])
    y_specific_heat.append(data1[19+20*(i),2])
    error_specific_heat.append(data1[19+20*(i),3])
    y_chi.append(data2[19+20*(i),2])
    error_chi.append(data2[19+20*(i),3])
    y_M.append(data3[19+20*(i),2])
    error_M.append(data3[19+20*(i),3])

#plotto l'energia interna per particella teorica e quella ottenuta con la simulazione
plt.errorbar(x, y_internal_energy, yerr=error_internal_energy, fmt='-', label="Simulated gibbs U/N", color='green', ecolor='yellow', capsize=3)
plt.plot(T, e, label="Theoretical U/N", color='blue')
plt.title('Ising 1D, internal energy')
plt.xlabel('T')
plt.ylabel('U/N')
plt.legend()
plt.savefig('grafici/internal_energy')

plt.clf()

#cv in funzione di T teorico
heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)
#plotto cv ottenuto con simulazione e cv teorico
plt.errorbar(x, y_specific_heat, yerr=error_specific_heat, fmt='-', label="Simulated gibbs cv", color='green', ecolor='yellow', capsize=3)
plt.plot(T, heat, label="Theoretical cv", color='blue')
plt.title('Ising 1D, Heat capacity')
plt.xlabel('T')
plt.ylabel('C')
plt.legend()
plt.savefig('grafici/heat_capacity')
plt.clf()

#suscettività in funzione di T teorica
X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
#plotto suscettività ottenuta con simulazione e teorica
plt.errorbar(x, y_chi, yerr=error_chi, fmt='-', label="Simulated gibbs chi", color='green', ecolor='yellow', capsize=3)
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

#plotto magnetizzazione ottenuta con simulazione e teorica
plt.errorbar(x, y_M, yerr=error_M, fmt='-', label="Simulated gibbs M", color='green', ecolor='yellow', capsize=3)
plt.plot(T, M, label="Theoretical M", color='blue')
plt.title('Ising 1D, magnetization M with h = 0.02')
plt.xlabel('T')
plt.ylabel('$M$')
plt.legend()
plt.savefig('grafici/magnetization')