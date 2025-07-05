/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//in questo codice si usa il metodo SA per trovare i parametri mu e sigma
//che minimizzano <H>
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  System SYS;
  SYS.initialize();//inizializzo il sistema
  SYS.initialize_properties();//inizializzo le proprietà
  SYS.block_reset(0);//inizializzo gli accumulatori di blocco
  //apro il file dove carico i valori di <H> e l'errore alla fine di ogni
  //simulazione MC
  ofstream coutf;
  coutf.open("../OUTPUT/SAevolution.dat");
  coutf << "# SAstep" << setw(20) << "<H>" << setw(20) << "error" << endl;
  //apro il file dove carico i valori di mu e sigma estratti e 
  //accettati durante la simulazione SA
  ofstream coutf1;
  coutf1.open("../OUTPUT/SAparameters.dat");
  coutf1 << "#mu" << setw(20) << "sigma" << endl;
  //variabili  di "deposito"
  double dep=0.0;
  double media=0.0;
  double sigma=0.0;
  double error_H=0.0;
  int SAstep=0;
  double mu_new=0.0;
  double sigma_new=0.0;
  //ciclo sulla temperatura. "raffreddamento"
  for(int u=0;u<100;u++){
  //ciclo di SA
  for(int y=0; y<SYS.getSAblocks(); y++){
    //ciclo MC
  for(int i=0; i < SYS.get_nbl(); i++){ //ciclo sui blocchi
    for(int j=0; j < SYS.get_nsteps(); j++){ //ciclo sugli steps di un blocco
      SYS.step();//effettuo uno step
      SYS.measure();//misura delle proprietà dopo lo step
    }
    SYS.averages(i+1);//calcolo medie di blocco, medie cumulative e incertezze di blocco
    SYS.block_reset(i+1);//resetto gli accumulatori di blocco
  }
  //se la mossa SA viene accettata, aggiorno le variabili <H>, errore di <H>, mu, sigma
  //la percentuale di mosse accettate dalla simulazione con Metropolis deve essere maggiore del 45%
  //secondo la "legge empirica del 50%"
    if(SYS.SAmove(dep,SYS.getH())==true&&SYS.get_acceptance()>0.45){
      dep=SYS.getH();
      media=SYS.get_media();
      sigma=SYS.get_sigma();
      error_H=SYS.get_error_H();
      SAstep++;
      coutf << SAstep << setw(20) << dep << setw(20) << error_H << endl;
      coutf1 << media << setw(20) << sigma << endl;
    };
    SYS.global_reset();//resetto gli accumulatori di medie cumulative e incertezze di blocco
    //propongo nuovi mu e sigma con controllo affinchè non siano minori o uguali a 0.2
    do{
      mu_new=media+SYS.Rannyu(-1.0,1.0)*0.05;
      sigma_new=sigma+SYS.Rannyu(-1.0,1.0)*0.05;
    }while(mu_new<=0.2||sigma_new<=0.2);
    SYS.set_media(mu_new);//propongo un nuovo mu
    SYS.set_sigma(sigma_new);//propongo un nuovo sigma
}
SYS.set_temperature(SYS.get_temperature()*0.9);//abbasso la temperatura
SYS.set_beta(SYS.get_temperature());//imposto beta di conseguenza
SYS.set_media(media);//imposto la nuova mu estratta
SYS.set_sigma(sigma);//imposto la nuova sigma estratta
}
  SYS.finalize();//finalizzo il sistema
  coutf.close();
  coutf1.close();
  //stampo a schermo il miglior valore di <h> e il suo errore, il miglior valore di mu e sigma e 
  //la temperatura finale della simulazione SA
  cout << dep << " +- " << error_H << endl;
  cout << media << " " << sigma << endl;
  cout << SYS.get_temperature() << endl;
  //carico su un file i migliori parametri mu e sigma ottenuti
   ofstream coutf2;
   coutf2.open("../OUTPUT/best_parameters.dat");
   coutf2 << "# mu" << setw(20) << "sigma" << endl;
   coutf2 << media << setw(20) << sigma << endl;
   coutf2.close();

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
