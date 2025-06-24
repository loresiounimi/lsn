#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "statistics.h"

using namespace std;

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

double Media (vector<double> v , int it){
  return (accumulate(v.begin(), v.end(), 0.0)/it);
  }

double DevStd (vector<double> v , int it){
  if(it==1)return 0;
  return sqrt((1.0/double(it-1))*((1.0/double(it))*sommaQuadrati(v)-pow(sommaElementi(v)/double(it),2)));
  }

double sommaQuadrati(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0,
        [](double accumulato, double elemento) {
            return accumulato + elemento * elemento;
        });
}

double sommaElementi(vector<double> v) {
    return accumulate(v.begin(), v.end(), 0.0);
}