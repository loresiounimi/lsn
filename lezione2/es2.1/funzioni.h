#ifndef __Funzione__
#define __Funzione__
#include <iostream>
#include <cmath>

using namespace std;
//classe madre
class funzione{
public:
  virtual double Eval(double x) const = 0;//metodo virtuale che implemento nelle classi figlie che restituisce f(x)
};

//classe figlia funzione g(x) per campionamento uniforme
class funzione1: public funzione{
public:
//implemento Eval
  double Eval(double x) const {
    return (M_PI/2)*cos((x*M_PI)/2);
  }
};
//classe figlia funzione g(x) per importance sampling ottenuta dividendo l'integranda per la p(x) da cui campiono le x_i
class funzione2: public funzione{
public:
//implemento Eval
  double Eval(double x) const {
    return (M_PI/2)*cos((x*M_PI)/2)/(-2*(x-1)); //
  }
};


#endif