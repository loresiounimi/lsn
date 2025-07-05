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
  const int _ndim = 3; // dimensionalità del sistema
  int _spin;            // spin della particella (+1 o -1)
  vec _x;               // posizione corrente
  vec _xold;            // posizione precedente (usata in moveback())
  vec _v;               // vettore delle velocità

public: // dichiarazione delle funzioni
  void initialize();                      // inizializza le proprietà della particella
  void translate(vec delta, vec side);   // sposta la particella entro il box di simulazione
  void flip();                           // flippa lo spin
  void moveback();                       // sposta la particella alla posizione precedente
  void acceptmove();                     // accetta la mossa proposta e muove la particella nella nuova posizione
  int  getspin();                        // restituisce lo spin della particella
  void setspin(int spin);                // imposta lo spin della particella
  double getposition(int dim, bool xnew);// restituisce la posizione della particella lungo una certa dimensione
  void   setposition(int dim, double position); // imposta la posizione della particella lungo una certa direzione
  void   setpositold(int dim, double position); // imposta la posizione precedente lungo una certa direzione
  double getvelocity(int dim);           // restituisce la velocità della particella lungo una certa direzione
  vec    getvelocity();                  // restituisce il vettore velocità della particella
  void   setvelocity(int dim, double velocity); // imposta la velocità di una particella lungo una certa dimensione
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
