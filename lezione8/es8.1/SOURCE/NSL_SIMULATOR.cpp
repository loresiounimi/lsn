/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//In questo esercizio ho aggiunto al file system un algoritmo per il calcolo di <H>.
//In questo codice ho sfruttato i migliori parametri ottenuti nell'esercizio 8.2 di simulated annealing per caricare su un file tutti i 
//valori estratti della posizione x campionati tramite metropolis secondo la distribuzione |psi|^2 per poter creare un grafico riempiendo
// un istogramma con tali valori e visualizzare la forma di |psi|^2
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  System SYS;
  SYS.initialize();//inizializzo il sistema
  SYS.initialize_properties();//inizializzo le proprietà
  SYS.block_reset(0);//inizializzo gli accumulatori di blocco
  //carico la miglior stima di mu e sigma ottenuta con simulated annealing
  double mu=0.0;
  double sigma=0.0;
  string header;
  ifstream parameters;
  //file dove si trova la miglior stima di mu e sigma
  parameters.open("../../es8.2/OUTPUT/best_parameters.dat");
  if (!parameters.is_open()) {
        cerr << "Errore nell'apertura del file best_parameters.dat" << endl;
        return 1;
    }
getline(parameters, header); 
parameters >> mu >> sigma; 
parameters.close();
SYS.set_media(mu);
SYS.set_sigma(sigma);
//apro un file su cui scrivere tutte le posizioni campionate da metropolis
//seguendo la distribuzione |psi|^2
ofstream coutf;
coutf.open("../OUTPUT/x_sampled.dat");
if (!coutf.is_open()) {
        cerr << "Errore nell'apertura del file x_sampled.dat" << endl;
        return 1;
    }
  for(int i=0; i < SYS.get_nbl(); i++){ //ciclo sui blocchi
    for(int j=0; j < SYS.get_nsteps(); j++){ //ciclo sugli steps di un blocco
      SYS.step();//effettuo uno step
      SYS.measure();//misura delle proprietà dopo lo step
      //carico le x campionate
      coutf << SYS.get_x(0) << endl;
        
    }
    SYS.averages(i+1);//calcolo medie di blocco, medie cumulative e incertezze di blocco
    SYS.block_reset(i+1);//resetto gli accumulatori di blocco
  }
  SYS.finalize();//finalizzo il sistema

  return 0;
}

//nota che il codice funziona anche con _npart>1

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
