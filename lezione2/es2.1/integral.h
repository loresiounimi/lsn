#ifndef __integral__
#define __integral__
#include <vector>
#include "funzioni.h"

 using namespace std;

double integraleMC (const funzione&, vector<double> , int);

double errore(vector<double> , int);
double sommaElementi(vector<double>);
#endif