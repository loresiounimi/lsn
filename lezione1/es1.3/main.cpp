#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <iomanip>
#include "random.h"
#include "statistics.h"

using namespace std;

int main (int argc, char *argv[]){
  //file su cui caricare i dati
  ofstream outputfile1("Buffon.txt");
  if(!outputfile1){
    cout << "Error: Unable to open file Buffon.txt" << endl;
    return 1;
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
    
  //calcolo il pigreco tramite la serie di Nilakantha
  double pi=0.0;
  for(int i=1;i<100000;i++){
      pi+=(4*pow(-1.0,double(i+1)))/(2.0*double(i)*(2.0*double(i)+1.0)*(2.0*double(i)+2.0));
    }
    pi+=3.0; 
    /*cout << setprecision(20) << pi << endl;
    cout << setprecision(20) << M_PI << endl;*/ //per confrontare quanto bene la serie converge a pigreco

    //esperimento di Buffon: ago di lunghezza L lanciato tra due righe orizzontali poste a distanza d
    int throws=10000; //lanci per blocco
    int N=100;//blocchi
    int hit=0;//accumulatore di quante volte l'ago colpisce le righe orizzontali della griglia
    double L=1.0;//lunghezza ago
    double d=1.2;//distanza tra due righe
    double center=0.0;//posizione dove cade il centro dell'ago, generata casualmente tra 0 e d
    vector <double> v;//accumulatore per medie di blocco per l stima di pigreco
    for(int j=0;j<N;j++){//ciclo sui blocchi
      hit=0.0;//inizializzo a 0 all'inizio di ogni blocco
      for(int i=0;i<throws;i++){//ciclo sui lanci
        center=rnd.Rannyu(0.0,d);//posizione del centro dell'ago
        //se l'ago interseca una delle due righe allora aggiorni hit
        if(center+(L/2.0)*sin(rnd.Rannyu(0.0,pi))>d||center-(L/2.0)*sin(rnd.Rannyu(0.0,pi))<0.0){
          hit++;
        }
        }
        v.push_back((2.0*L*throws)/(hit*d));//calcolo di pigreco per il singolo blocco
        //scrivo su file il numero di lanci, le medie cumulative e l'incertezza ottenuta con data blocking
        outputfile1 << (j+1)*throws << " " << sommaElementi(v)/(j+1) << " " << DevStd(v,j+1) << endl;
      }
      

    rnd.SaveSeed();
  return 0;
}