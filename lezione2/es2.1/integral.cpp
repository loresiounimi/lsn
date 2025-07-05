#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "funzioni.h"
#include "integral.h"

using namespace std;
//calcolo dell'integrale con metodo montecarlo 
double integraleMC (const funzione& f, vector<double> v, int n){
  double sum=0.0;
  for (int i=0;i<n;i++){
    sum+=f.Eval(v[i]);
  }
  return sum/n;
}
//calcolo dell'errore con l'importance sampling
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
//funzione che restituisce la somma degli elementi di un vector
double sommaElementi(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0);
}