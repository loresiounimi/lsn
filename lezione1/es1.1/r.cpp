#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include "statistics.h"
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
  //creo i file su cui caricare i dati
  ofstream outputfile1("<r>.txt");
  ofstream outputfile2("r_sigma.txt");
  if(!outputfile1){
    cout << "Error: Unable to open file <r>.txt" << endl;
    return 1;
    } 
  if(!outputfile2){
    cout << "Error: Unable to open file r_sigma.txt" << endl;
    return 1;
    } 
  
  int nmax= 100000; //numero massimo di MC steps
  int N=100; //numero di blocchi MC
  int l; 
  l=nmax/N; //numero steps MC per blocco
  //inizializzo il seed
  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.out");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  vector <vector <double>> v(2);//contenitore per memorizzare medie di blocco
  int it=0;//variabile che tiene conto dei MC steps fatti
  //variabili "deposito" per memorizzare i dati delle l simulazioni di un blocco
  double dep = 0.0;
  double dep2 = 0.0;
  for(int j=0; j<N; j++){
    //inizializzo i dep a zero per avere medie di blocco indipendenti
    dep = 0.0;
    dep2= 0.0;
    for (int i=0;i<l;i++){
      dep+=rnd.Rannyu();
      dep2+=pow(rnd.Rannyu()-0.5,2);
      it++;
    }
    v[0].push_back(dep/double(l));  //media di blocco di <r>
    v[1].push_back(dep2/double(l));  //media di blocco di sigma
    //carico i dati sui file
    outputfile1 << sommaElementi(v[0])/(j+1) << " " << it << " " << DevStd(v[0],j+1) << endl;
    outputfile2 << sommaElementi(v[1])/(j+1) << " " << it << " " << DevStd(v[1],j+1) << endl;
    // Output su file: media cumulativa | passi Monte Carlo | errore statistico sul blocco
  }
  rnd.SaveSeed();
  return 0;
} 