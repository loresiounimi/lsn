/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//la classe particle descrive un oggetto particella con le sue posizioni attuali x,y,z, le sue posizioni al tempo precedente x,y,z, le sue velocità vx, vy, vz e il suo spin
#ifndef __Particle__
#define __Particle__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Particle {

private:
  const int _ndim = 1; // dimensionalità del sistema
  vec _x;               // posizione corrente
  vec _xold;            // posizione precedente (usata in moveback())

public: // dichiarazione delle funzioni
  void initialize();                      // inizializza le proprietà della particella
  void translate(vec delta);              // sposta la particella 
  void moveback();                       // sposta la particella alla posizione precedente
  void acceptmove();                     // accetta la mossa proposta e muove la particella nella nuova posizione
  double getposition(int dim, bool xnew);// restituisce la posizione della particella lungo una certa dimensione
  void   setposition(int dim, double position); // imposta la posizione della particella lungo una certa direzione
  void   setpositold(int dim, double position); // imposta la posizione precedente lungo una certa direzione
  double pbc(double position, double side);     // applica le pbc

};

#endif // __Particle__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
