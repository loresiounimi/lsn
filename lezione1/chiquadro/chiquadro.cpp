#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"
#include "statistics.h"

using namespace std;

int main (int argc, char *argv[]){
  //apro i file su cui scrivere i dati
  ofstream outputfile1("chiquadro_singoli.txt");
  if(!outputfile1){
    cout << "Error: Unable to open file chiquadro_singoli.txt" << endl;
    return 1;
    } 
  ofstream outputfile2("chiquadro.txt");
  if(!outputfile1){
    cout << "Error: Unable to open file chiquadro.txt" << endl;
    return 1;
    } 
  
  int M= 100;//numero di sottointervalli in cui è diviso l'intervallo [0,1]
  int N=10000;//numero di estrazioni di numeri casuali con rnd
  int l=100;//numero di chiquadro da calcolare
  double dep=0;
  double step=0.01;//larghezza dei sottointervalli

  //inizializzo il generatore di numeri casuali
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
  vector <int> x(M);//vettore che accumula le occorrenze dei sottointervalli
  vector <double> chiquadro;//vettore che accumula i chiquadro ottenuti dai 100 diversi cicli
  for(int j=0; j<l; j++){//ciclo sui 100 chiquadro
    for(int i=0; i<N; i++){//ciclo sul sungolo chiquadro
      dep=rnd.Rannyu();
      for(int h=0;h<M;h++){//riempimento del vettore delle occorrenze
      if(dep>(h*step) && dep<((h+1)*step)){
        x[h]++;
      }
      }
    }
    chiquadro.push_back(ChiQuadro(x, M, double(N/M) ));//calcolo del j-esimo chiquadro
    //scrittura su file
    outputfile1 << chiquadro[j] << endl;//chi quadro singolo
    outputfile2 << sommaElementi(chiquadro)/(j+1) << " " << DevStd(chiquadro,j+1) << " " << j+1 << endl;//medie cumulative e incertezza con data blocking
    fill(x.begin(), x.end(), 0);//x è il vettore che memorizza e accumula le occorrenze, alla fine di ogni ciclo sul singolo chiquadro lo inizializzo a 0
  }
  rnd.SaveSeed();//salvo il seed
  return 0;
}