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
  Random rnd;
  int M= 100;
  int N=10000;
  double dep=0;
  double step=0.01;
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
  vector <int> x(M);
  vector <double> chiquadro;
  for(int j=0; j<100; j++){
    for(int i=0; i<N; i++){
      dep=rnd.Rannyu();
      for(int h=0;h<M;h++){
      if(dep>(h*step) && dep<((h+1)*step)){
        x[h]++;
      }
      }
    }
    chiquadro.push_back(ChiQuadro(x, M, double(N/M) ));
    outputfile1 << chiquadro[j] << endl;
    outputfile2 << sommaElementi(chiquadro)/(j+1) << " " << DevStd(chiquadro,j+1) << " " << j+1 << endl;
    fill(x.begin(), x.end(), 0);
  }
  rnd.SaveSeed();
  return 0;
}