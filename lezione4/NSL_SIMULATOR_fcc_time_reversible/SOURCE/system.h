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

//system è una classe in grado di eseguire vari tipi di simulazione molecolare e misurare alcune proprietà medie della simulazione nella sua evoluzione temporale
//in questo esercizio è effettuatta una simulazione di MD in modello LJ

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> // exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 3;  // dimensionalità del sistema
  bool _restart;        // flag che indica se la simulazione è un restart
  int _sim_type;        // tipo di simulazione (es. lennard-jones, ising)
  int _npart;           // numero di particelle
  int _nblocks;         // numero di blocchi della simulazione
  int _nsteps;          // numero di passi della simulazione per blocco
  int _nattempts;       // numero di mosse tentate
  int _naccepted;       // numero di mosse accettate
  double _temp, _beta;  // temperatura e beta (1/kbT)
  double _rho, _volume; // densità e volume del sistema
  double _r_cut;        // raggio di cutoff per interazioni tra particelle
  double _delta;        // ampiezza massima del passo di spostamento
  double _J, _H;        // parametri dell'hamiltoniana di ising
  vec _side;            // dimensioni della scatola
  vec _halfside;        // metà delle dimensioni della scatola
  Random _rnd;          // generatore di numeri casuali
  field<Particle> _particle; // vettore di oggetti particella che rappresentano il sistema
  vec _fx, _fy, _fz;    // forze sulle particelle nelle direzioni x, y e z
  
  // proprietà
  int _nprop; // numero di proprietà misurate
  bool _measure_penergy, _measure_kenergy, _measure_tenergy; // flag per misurare le varie energie
  bool _measure_temp, _measure_pressure, _measure_gofr;      // flag per misurare temperatura, pressione e funzione di distribuzione radiale
  bool _measure_magnet, _measure_cv, _measure_chi;           // flag per misurare magnetizzazione, capacità termica e suscettività
  bool _measure_pofv;                                        // flag per misurare la distribuzione del modulo della velocità
  int _index_penergy, _index_kenergy, _index_tenergy;        // indici per accedere alle proprietà relative all’energia potenziale, cinetica e totale
  int _index_temp, _index_pressure, _index_gofr;             // indici per accedere a temperatura, pressione e g(r)
  int _index_magnet, _index_cv, _index_chi;                  // indici per accedere a magnetizzazione, capacità termica e suscettività
  int _index_pofv;                                           // indice per accedere alla distribuzione del modulo della velocità
  int _n_bins;           // numero di bin per la funzione di distribuzione radiale
  int _n_bins_v;         // numero di bin per la distribuzione del modulo della velocità
  double _bin_size;      // ampiezza dei bin per g(r)
  double _bin_size_v;    // ampiezza dei bin per la distribuzione del modulo della velocità
  double _vtail, _ptail; // correzioni di coda per energia e pressione
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
  void write_XYZ(int nconf);  // scrive la configurazione attuale in formato xyz durante la simulazione
  void read_configuration();  // legge la configurazione del sistema da file
  void initialize_velocities();// inizializza le velocità delle particelle
  void step();                // esegue un passo di simulazione
  void block_reset(int blk);  // resetta le medie del blocco corrente
  void measure();             // misura le proprietà del sistema
  void averages(int blk);     // calcola le medie delle proprietà nel blocco corrente
  double error(double acc, double acc2, int blk); // calcola l’errore statistico
  void move(int part);        // muove una particella
  bool metro(int part);       // esegue lo step di metropolis per una particella
  double pbc(double position, int i); // applica le condizioni periodiche alle coordinate
  int pbc(int i);             // applica le condizioni periodiche ai siti di spin
  void Verlet();              // esegue l’integrazione con l’algoritmo di verlet
  double Force(int i, int dim); // calcola la forza sulla particella i nella direzione dim
  double Boltzmann(int i, bool xnew); // calcola il fattore di boltzmann per l’accettazione metropolis
  //funzioni aggiunte per questo esercizio di reversibilità temporale
  void read_configuration_timereversible();  // legge la configurazione del sistema da file
  void finalize_timereversible();            // finalizza il sistema dopo la simulazione in avanti nel tempo
  void write_configuration_timereversible();  // scrive la configurazione finale su dei file preparando la simulazione indietro nel tempo

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
