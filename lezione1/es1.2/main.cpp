#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
  //file su cui scriverò i dati
  string isto[3]= {"standardisto.txt", "expisto.txt", "cauchyisto.txt"};
  ofstream istowrite[3];
  for(int i=0;i<3;i++){
    istowrite[i].open(isto[i]);
  if(!istowrite[i]){
    cout << "Error: Unable to open the file" << endl;
    return 1;
    } 
    } 
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
  int nmax= 10000;//realizzazioni di S_N
  double dep=0.0;
  vector <int> s{1,2,10,100};//valore di N in S_N

  //vector di dimensione 3 (1 per ogni distribuzione da analizzare) che a loro volta contengono 4 vector(1 per ogni valore di N) che a loro volta contengono vector di double che svolgono la funzione di accumulatori di occorrenze
  //per gli istogrammi da realizzare
  vector <vector<vector<double>>> v(3);

  for(int y=0; y<3; y++){//ciclo sulle distribuzioni
    v[y].resize(s.size());//assegno la dimensione dei primi vector contenuti nel vector v
  for(int j=0; j<(int)s.size(); j++){//ciclo sui 4 possibili valori di N
  for(int i=0; i<nmax; i++){//ciclo sulle 10000 realizzazioni di S_N
    for(int h=0; h<s[j]; h++){
      if(y==0){
        dep+=rnd.Rannyu();//distribuzione uniforme
        }else if(y==1){
        dep+=rnd.exp(1.0);//distribuzione esponenziale
        }else if(y==2){
        dep+=rnd.Cauchy(0.0,1.0);//distribuzione di cauchy-lorentz
        }else{
        break;
        }
      }
    v[y][j].push_back(dep/s[j]);//memorizzo S_N
    dep=0.0;
    }
    }
  //scrivo sui file
  for(int i=0; i<nmax; i++){
    istowrite[y] << v[y][0][i] << " " << v[y][1][i] << " " << v[y][2][i] << " " << v[y][3][i] << endl;
    }
    }
    rnd.SaveSeed();
  return 0;
}