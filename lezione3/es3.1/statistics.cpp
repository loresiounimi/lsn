#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "funzioni.h"
#include "statistics.h"

using namespace std;

//funzione per il calcolo dell'incertezza con il data blocking
double errore(vector<double> v , int n){
  if(n<=1){return 0;}else{
  double sum1=(accumulate(v.begin(), v.end(), 0.0))/n;
  double sum2=(accumulate(v.begin(), v.end(), 0.0, 
  [](double acc, double val) {
      return acc + val * val;
  }
  //con accumulate sommo gli elementi del vettore v partendo dal valore iniziale della somma pari a 0.0
  //implemento per sum2 una lambda function per poter sommare gli elementi di v al quadrato
))/n;
return sqrt((sum2-pow(sum1,2))/(n-1));
}
}
//funzione che somma gli elementi di un vettore
double sommaElementi(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0);
}