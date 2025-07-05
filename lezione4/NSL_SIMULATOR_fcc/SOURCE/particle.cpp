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
   _spin = 1;
   _x.resize(_ndim);
   _xold.resize(_ndim);
   _v.resize(_ndim);
   return;
}
//sposta la particella di una quantità delta lungo ogni direzione e rispettando le pbc
void Particle :: translate(vec delta, vec side){
   for(unsigned int i=0; i<_ndim; i++){
     _x(i) = pbc(_x(i) + delta(i), side(i));
   }
}
//flippa lo spin della particella
void Particle :: flip(){
   _spin = -1*this->getspin();
}
//fa tornare la particella alla posizione precdente
void Particle :: moveback(){
   _x = _xold;
}
//accetta la mossa proposta di spostamento della posizione della particella
void Particle :: acceptmove(){
   _xold = _x;
}
//restitusce il valore dello spin della particella
int Particle :: getspin(){
   return _spin;
}
//setta lo spin della particella
void Particle :: setspin(int spin){
   _spin = spin;
   return;
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
//restituisce una delle 3 velocità spaziali
double Particle :: getvelocity(int dim){
   return _v(dim);
}
//restituisce il vettore delle 3 velocità spaziali
vec Particle :: getvelocity(){
   return _v;
}
//setta una delle 3 velocità spaziali
void Particle :: setvelocity(int dim, double velocity){
   _v(dim) = velocity;
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
