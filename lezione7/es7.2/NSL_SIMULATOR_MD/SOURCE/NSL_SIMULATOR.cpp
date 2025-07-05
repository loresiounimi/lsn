/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//simulazione MD con potenziale Lennard Jones
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();//inizializzo il sistema
  SYS.initialize_properties();//inizializzo le proprietà
  SYS.block_reset(0);//inizializzo gli accumulatori di blocco

  for(int j=0; j < 5000; j++){ //ciclo per equilibrare il sistema
    SYS.step();
    }
  for(int i=0; i < SYS.get_nbl(); i++){ //ciclo sui blocchi
    for(int j=0; j < SYS.get_nsteps(); j++){ //ciclo sugli steps di un blocco
      SYS.step();//effettuo uno step
      SYS.measure();//misura delle proprietà dopo lo step
    }
    SYS.averages(i+1);//calcolo medie di blocco, medie cumulative e incertezze di blocco
    SYS.block_reset(i+1);//resetto gli accumulatori di blocco
  }
  SYS.finalize();//finalizzo il sistema

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
