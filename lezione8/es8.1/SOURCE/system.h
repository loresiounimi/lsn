/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 1;  // dimensionalità del sistema
  bool _restart;        // flag che indica se la simulazione è un restart
  int _npart;           // numero di particelle
  int _nblocks;         // numero di blocchi della simulazione
  int _nsteps;          // numero di passi della simulazione per blocco
  int _nattempts;       // numero di mosse tentate
  int _naccepted;       // numero di mosse accettate
  double _temp, _beta;  // temperatura e beta (1/kbT)
  double _delta;        // ampiezza massima del passo di spostamento
  double _media, _sigma; //parametri di psi
  Random _rnd;          // generatore di numeri casuali
  field<Particle> _particle; // vettore di oggetti particella che rappresentano il sistema
  
  // Propietà
  int _nprop; // numero di proprietà misurate
  bool _measure_H;                                           //flag per misurare il valor medio di H
  int _index_H;           //indice per accedere alle proprietà di H
  vec _block_av;         // medie a blocchi delle proprietà
  vec _global_av;        // medie globali delle proprietà
  vec _global_av2;       // quadrati delle medie globali delle proprietà
  vec _average;          // medie finali delle proprietà
  vec _measurement;      // misure correnti delle proprietà

public: // dichiarazione delle funzioni
  int get_nbl();              // restituisce il numero di blocchi
  int get_nsteps();           // restituisce il numero di passi per blocco
  void initialize();          // inizializza le proprietà del sistema
  void initialize_properties();// inizializza le proprietà da misurare
  void finalize();            // finalizza la simulazione e libera eventuali risorse
  void write_configuration(); // scrive la configurazione finale del sistema in formato xyz
  void read_configuration();  // legge la configurazione del sistema da file
  void step();                // esegue un passo di simulazione
  void block_reset(int blk);  // resetta le medie del blocco corrente
  void measure();             // misura le proprietà del sistema
  void averages(int blk);     // calcola le medie delle proprietà nel blocco corrente
  double error(double acc, double acc2, int blk); // calcola l’errore statistico
  void move(int part);        // muove una particella
  bool metro(int part);       // esegue lo step di metropolis per una particella
  //funzioni aggiunte per questo esercizio 
  void set_temperature(double t); // funzione che imposta la temperatura
  double get_temperature(); //funzione che restituisce la temperatura
  void set_beta(double t); //funzione che imposta beta
  void global_reset(); //funzione che azzera tutti gli accumulatori
  void set_media(double media); //funzione che setta il valore di mu
  double get_media(); //funzione che restituisce il valore di mu
  void set_sigma(double sigma); //funzione che setta il valore di sigma
  double get_sigma(); //funzione che restituisce il valore di sigma
  double get_x(int i); //funzione che restituisce la posizione 1D (per convenzione x) della particella i-esima
  double PSI(int i, bool xnew);  //funzione che calcola il valore di PSI(x)
  double PSI_seconda(int i, bool xnew); //funzione che calcola il modulo quadro di PSI(x)
};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
