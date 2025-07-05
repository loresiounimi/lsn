/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//simulazione di un modello di Ising 1D con algoritmo di Gibbs
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();//inizializzo il sistema
  SYS.initialize_properties();//inizializzo le proprietà
  for(int l=0;l<26;l++){//ciclo sulle temperature
    
    SYS.block_reset(0);//inizializzo gli accumulatori di blocco a zero all'inizio di ogni simulazione per ogni T
  for(int j=0; j < 10000; j++){ //ciclo per far equilibrare il sistema in modo da non avere transitori
    SYS.step();
    }
  for(int i=0; i < SYS.get_nbl(); i++){ //ciclo sui blocchi
    for(int j=0; j < SYS.get_nsteps(); j++){ //ciclo sugli steps di un blocco
      SYS.step();//effettuo uno step
      SYS.measure();//misura delle proprietà dopo lo step
      
    }
    SYS.averages(i+1);//calcolo medie di blocco, cumulative e incertezze con data blocking
    SYS.block_reset(i+1);//resetto gli accumulatori di blocco
  }
  SYS.set_temperature(SYS.get_temperature()+0.1);//imposto la nuova temperatura
  SYS.set_beta(SYS.get_temperature());//imposto il nuovo beta
  SYS.global_reset(); //resetto gli accumulatori delle medie cumulative e delle incertezze per iniziare una simulazione con una nuova T
}
  SYS.finalize();

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
