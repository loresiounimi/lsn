#ifndef __Funzione__
#define __Funzione__
#include <iostream>
#include <cmath>
#include "random.h"
using namespace std;

//classe madre
class funzione{
public:
  virtual double Eval(double x) const=0; //metodo virtuale che implemento nelle classi figlie per restituire il valore f(x)
  
};

//Implemento la funzione di Wiener che descrive un GBM. E' la classe figlia
class wiener : public funzione {
  private:
      double sigma;     // Volatilità
      double interest;  // Tasso di interesse
      double S0;        // prezzo iniziale
      double tprec;    // tempo al passo precedente
      Random& rnd;     // Generatore di numeri casuali (passato per riferimento, cioè uso l'oggetto originale)
  
  public:
      // Costruttore per inizializzare i parametri
      wiener(double sigma, double interest, double S0, Random& rnd, double tprec)
          : sigma(sigma), interest(interest), S0(S0), rnd(rnd), tprec(tprec) {}
      //funzione per settare S0
      void setS0(double S0new){
        S0=S0new;
      }
      //funzione per settare t al passo precedente
      void settprec(double tprecnew){
        tprec=tprecnew;
      }
      // Implementazione di Eval
      double Eval(double x) const override{
          return S0 * exp((interest - pow(sigma, 2) / 2) * (x-tprec) + sigma * rnd.Gauss(0.0, 1.0) * sqrt(x-tprec));
      }
  };

#endif