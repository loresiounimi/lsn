#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}
// Costruttore di default, non esegue alcuna azione

Random :: ~Random(){}
// Distruttore di default, non esegue alcuna azione

void Random :: SaveSeed(){
   // Questa funzione salva lo stato attuale del generatore di numeri casuali nel file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED\t" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
   } else cerr << "PROBLEMA: impossibile aprire random.out" << endl;
   WriteSeed.close();
   return;
}

double Random :: Gauss(double mean, double sigma) {
   // Questa funzione genera un numero casuale con distribuzione gaussiana di media `mean` e deviazione standard `sigma`
   double s = Rannyu();
   double t = Rannyu();
   double x = sqrt(-2.*log(1.-s)) * cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // Questa funzione genera un numero casuale nell'intervallo [min, max)
   return min + (max - min) * Rannyu();
}

double Random :: Rannyu(void){
  // Questa funzione genera un numero casuale nell'intervallo [0,1)
  const double twom12 = 0.000244140625;  // 2^-12
  int i1, i2, i3, i4;
  double r;

  // Calcoli per aggiornare lo stato interno del generatore
  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;

  // Aggiornamento dei parametri
  l4 = i4 % 4096;
  i3 = i3 + i4 / 4096;
  l3 = i3 % 4096;
  i2 = i2 + i3 / 4096;
  l2 = i2 % 4096;
  l1 = (i1 + i2 / 4096) % 4096;

  // Calcolo del numero casuale tra 0 e 1
  r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  // Questa funzione imposta il seme e i parametri del generatore di numeri casuali
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random :: exp(double lambda){
  //questa funzione genera un numero casuale secondo una distribuzione esponenziale con parametro lambda
  return -log(1.0-Rannyu())/lambda;
}

double Random :: Cauchy(double mu, double gamma){
  //questa funzione genera un numero casuale secondo la distribuzione di Cauchy-Lorentz con parametri mu e gamma
  return (mu+gamma*tan((M_PI)*(Rannyu()-0.5)));
}