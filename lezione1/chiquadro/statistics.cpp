#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "statistics.h"

using namespace std;
//funzione che calcola il chiquadro
double ChiQuadro (vector<int> v , int it, double m){
    double dep=0.0;
    if(int(v.size())!=it){
        cout << "errore" << endl;
        return 0;
        }
    for(int i=0;i<it;i++){
        dep+=pow(double(v[i]-m),2);
        }
    return dep/m;
  }
//funzione per calcolare la media dei dati contenuti in un vettore di dimensione it
double Media (vector<double> v , int it){
  return (accumulate(v.begin(), v.end(), 0.0)/it);
  }
//funzione per calcolare l'errore MC su dei dati contenuti nel vettore in funzione del numero di blocchi it della simulazione MC 
double DevStd (vector<double> v , int it){
  if(it==1)return 0;
  return sqrt((1.0/double(it-1))*((1.0/double(it))*sommaQuadrati(v)-pow(sommaElementi(v)/double(it),2)));
  }

  //funzione che somma i valori di un vettore al quadrato
double sommaQuadrati(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0,
        [](double accumulato, double elemento) {
            return accumulato + elemento * elemento;
        });
         //con accumulate sommo gli elementi del vettore v partendo dal valore iniziale della somma pari a 0.0
         //implemento per sum2 una lambda function per poter sommare gli elementi di v al quadrato
}
//funzione che somma i valori di un vettore
double sommaElementi(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0);
}