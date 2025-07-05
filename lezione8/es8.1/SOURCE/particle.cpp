/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <math.h>
#include "particle.h"

using namespace std;

//inizializza la particella con la dimensionalità del sistema e lo spin
void Particle :: initialize(){

   _x.resize(_ndim);
   _xold.resize(_ndim);
   return;
}
//sposta la particella di una quantità lungo ogni direzione. In questa simulazione
//non si applicano le pbc
void Particle :: translate(vec delta){
   for(unsigned int i=0; i<_ndim; i++){
     _x(i) = _x(i) + delta(i);
   }
}

//fa tornare la particella alla posizione precdente
void Particle :: moveback(){
   _x = _xold;
}

//accetta la mossa proposta di spostamento della posizione della particella
void Particle :: acceptmove(){
   _xold = _x;
}

//restituisce una delle tre coordinate spaziali
double Particle :: getposition(int dim, bool xnew){
   if(xnew) return _x(dim);
   else return _xold(dim);
}
//setta una delle 3 coordinate spaziali
void Particle :: setposition(int dim, double position){
   _x(dim) = position;
   return;
}

//setta una delle 3 coordinate spaziali precdenti
void Particle :: setpositold(int dim, double position){
   _xold(dim) = position;
   return;
}
//applica le pbc
double Particle :: pbc(double position, double side){
  return position - side * rint(position / side);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
