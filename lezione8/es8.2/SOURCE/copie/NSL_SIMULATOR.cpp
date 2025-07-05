/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("../INPUT/Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("../OUTPUT/seed.out");
  if (input.is_open()){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
  } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  input.close();
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  /*for(int j=0; j < 10; j++){ //loop over steps to equilibrate the system
    SYS.step();
    }*/
   ofstream coutf;
   coutf.open("../OUTPUT/SAevolution.dat");
   coutf << "# SAstep" << setw(20) << "<H>" << setw(20) << "error" << endl;
   ofstream coutf1;
   coutf1.open("../OUTPUT/SAparameters.dat");
   coutf1 << "# SAstep" << setw(20) << "mu" << setw(20) << "sigma" << endl;
   double dep=0.0;
   double media=0.0;
   double sigma=0.0;
   double error_H=0.0;
   int SAstep=0;
   for(int u=0;u<200;u++){
    //cout << SYS.get_temperature() << endl;
  for(int y=0; y<SYS.getSAblocks(); y++){
    
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    //cout << SYS.get_media() << endl;
    //cout << SYS.getH() << " H" << endl;
    if(SYS.SAmove(dep,SYS.getH())==true&&SYS.get_acceptance()>0.45){
      dep=SYS.getH();
      media=SYS.get_media();
      sigma=SYS.get_sigma();
      error_H=SYS.get_error_H();
      SAstep++;
      coutf << SAstep << setw(20) << dep << setw(20) << error_H << endl;
      coutf1 << SAstep << setw(20) << media << setw(20) << sigma << endl;
    };
    //cout << media << " " << dep << " " << error_H << endl;
    SYS.global_reset();
    SYS.set_media(fabs(rnd.Rannyu(0.2, 1.2)));
    //cout << SYS.get_media() << endl;
    SYS.set_sigma(fabs(rnd.Rannyu(0.2,1.2)));
    //cout << SYS.get_media() << " random" << endl;
}
SYS.set_temperature(SYS.get_temperature()*0.9);
SYS.set_beta(SYS.get_temperature());
SYS.set_media(media);
SYS.set_sigma(sigma);
}
  //cout << SYS.getH() << endl;
  SYS.finalize();
  coutf.close();
  coutf1.close();
  cout << dep << " +- " << error_H << endl;
  cout << media << " " << sigma << endl;
  cout << SYS.get_temperature() << endl;
  rnd.SaveSeed();
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
