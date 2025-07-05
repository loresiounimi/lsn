/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // esegue un passo della simulazione
  if(_sim_type == 0) this->Verlet();  // esegue un passo di simulazione MD
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // esegue un passo di simulazione MC su una particella scelta a caso
  _nattempts += _npart; //aggiorna il numero di tentativi eseguiti sul sistema
  return;
}

void System :: Verlet(){ //esegue l'algoritmo di Verlet per aggiornare la posizione delle particelle del sistema 
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //forza agente sulla particella i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //algoritmo di Verlet
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;//aggiorna il numero di mosse accettate. In una simulazione MD con Verlet esso è pari al numero di tentativi eseguiti
  return;
}

double System :: Force(int i, int dim){ //forza agente sulla particella i lungo la direzione dim. E' calcolata attraverso il potenziale di L-J
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

//questa funzione non verrà usata in questo esercizio
void System :: move(int i){ // Propone una mossa per simulazione MC
  if(_sim_type == 3){ //Gibbs sampler for Ising
    // TO BE FIXED IN EXERCISE 6
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}
//questa funzione non verrà usata in questo esercizio
bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}
//questa funzione non verrà usata in questo esercizio
double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // applica le pbc
  return position - _side(i) * rint(position / _side(i));
}
//questa funzione non verrà usata in questo esercizio
int System :: pbc(int i){ // applica le pbc per gli spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(){ // inizializza il sistema in base al contenuto dei file di input che si trovano nella cartella INPUT

  int p1, p2; // inizializza il rnd
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; 
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Setta la prima riga del file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // inizia a leggere il contenuto del file ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");//apre il file di OUTPUT su cui scrivo poi i dati di input
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){ //inizializza il sistema in base al tipo di simulazione
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl; //simulazione di MD con potenziale LJ, quella che èutilizzata nell'esercizio corrente
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){ //restart dice al sistema se generare valori per le velocità e dunque per le posizioni precdenti oppure se leggere le posizioni precedenti da un file già esistente
      input >> _restart;
    } else if( property == "TEMP" ){//memorizza la temperatura e il suo inverso beta
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){//memorizza il numero di particelle del sistema e inizializza con questo valore la dimensione dei vettori delle forze e il vettore particella
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ //inizializza tutte le NPART particelle
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // setta randomicamente lo spin delle particelle, non usato in questo esercizio
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){//inizializza la densità e pertanto il volume 
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim); 
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){//inizializza rcut, cioè la distanza entro cui una particella risente dell'interazione con un'altra particella
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){//passo della simulazione MD
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){//numero di blocchi della simulazione
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){//numero di steps per blocco
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;//scrive sul file di output che la letture dei dati di input è andata abuon fine
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();//legge le posizioni iniziali da un file già esistente
  if(_sim_type==0) this->initialize_velocities();//inizializzo le velocità iniziali
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_velocities(){//inizializzo le velocità
  double xold, yold, zold;
  if(_restart){//Verlet per calcolare le velocità delle particelle utilizza la posizione a t+deltat e a t-deltat, dunque se esiste un file con le config delle particelle a t-delta t e cioè a t precedente, le legge
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/conf-1.xyz");//apre il file da cui leggere le posizioni precedenti e le memorizza con le pbc
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){//se il numero di particelle del sistema non corrisponde al numero di particelle nel file il codice dà errore
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      for(int i=0; i<_npart; i++){//setta xold secondo le pbc
        cinf >> particle >> xold >> yold >> zold; 
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  } else {//se non è presente il file con le xold, le calcola assegnando casualmente le velocità secondo una distribuzione di boltzmann e dalle velocità assegna le xold 
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    //assegna casualmente le velocità secondo la distribuzione di maxwell boltzmann alla temperatura temp e le somma per calcolare la velocità media del CM del sistema
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);//calcola la velocità del CM del sistema
    double sumv2 = 0.0, scalef;
    //sottrae la velocità del CM del sistema: questo annulla il momento totale del sistema, come richiesto dalla fisica statistica per sistemi isolati. 
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // calcola scalef
    //assegna le velocità in modo che l’energia cinetica media corrisponda esattamente alla temperatura desiderata ((3/2)kT, con k=1 in unità ridotte).
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
    //calcola e assegna xold in base alle velocità iniziali
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System :: initialize_properties(){ // inizializza i data membri utilizzati per la misura delle proprietà del sistema
  //ogni proprietà da misurare ha il suo indice, serve per poter accedere ai vettori che svolgono le operazioni di media e data blocking
  string property;
  int index_property = 0;
  _nprop = 0;//inizialmente il numero di proprietà è pari a zero

  _measure_penergy  = false; //definirò nel file properties le proprietà da misurare
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_pofv     = false;
  ifstream input("../INPUT/properties.dat");//apro il file contenente le proprietà da misurare e inizio a leggerlo
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){//energia potenziale
        ofstream coutp("../OUTPUT/potential_energy.dat");//apro il file su cui scriverò l'energia potenziale nella simulazione MD
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 0.0; // non serve in questo esercizio
      } else if( property == "KINETIC_ENERGY" ){//energia cinetica
        ofstream coutk("../OUTPUT/kinetic_energy.dat");//apro il file su cui scriverò l'energia cinetica nella simulazione MD
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
       
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){//energia totale
        ofstream coutt("../OUTPUT/total_energy.dat");//apro il file su cui scriverò l'energia totale  nella simulazione MD
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){//temperatura
        ofstream coutte("../OUTPUT/temperature.dat");//apro il file su cui scriverò la temperatura nella simulazione MD
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){//pressione
        ofstream coutpr("../OUTPUT/pressure.dat");//apro il file su cui scriverò la pressione nella simulazione MD
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 0.0; // non serve in questo esercizio
      } else if( property == "GOFR" ){//non serve in questo esercizio
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "POFV" ){//distribuzione delle velocità del sistema. le suddivido in bins 
        if(_sim_type > 0){
          cerr << "PROBLEM: DOES NOT MAKE SENSE COMPUTING POFV FOR MC" << endl;//si calcola solo per simulazioni MD 
          exit(EXIT_FAILURE);
        }
        ofstream coutpv("../OUTPUT/pofv.dat");//apro il file su cui scriverò i bins delle velocità nella simulazione MD
        coutpv << "# VELOCITY:     ACTUAL_POFV:     AVE_POFV:        ERROR:" << endl;
        coutpv.close();
        input>>_n_bins_v;//carico il numero di bins
        _nprop += _n_bins_v;
        _bin_size_v = 4.0*sqrt(_temp)/(double)_n_bins_v; // setto la dimensione dei bin in modo intelligente
        _measure_pofv = true;
        _index_pofv = index_property;
        index_property += _n_bins_v;
      } else if( property == "MAGNETIZATION" ){//non serve in questo esercizio
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){//non serve in questo esercizio
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){//non serve in questo esercizio
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){//fine lettura del file properties
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  //in base al numero di proprietà da misurare, assegno la dimensione dei vettori _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System :: finalize(){//finalizzo il sistema
  this->write_configuration();//scrivo le posizioni finali delle particelle
  _rnd.SaveSeed();//salvo il seed
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);//scrivo sul file di output che la simulazione è finita
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// scrivo la posizione finale delle particelle in un file.xyz nella cartella ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){//se non è una simulazione di ISING
    coutf.open("../OUTPUT/CONFIG/config.xyz");//apro il file
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){//carico le posizioni finali
        coutf << "LJ" << "  " 
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");//carico le posizioni precedenti alle posizioni finali
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {//non serve in questo esercizio
    coutf.open("../OUTPUT/CONFIG/config.spin");//carico gli spins finali delle particelle
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// scrive configurazione su un file xyz, serve per visualizzare l'evoluzione del sistema in 3D. non usato in questo codice
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

// leggo le posizioni iniziali da un file in ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config.xyz");//apro il file
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){//carico le posizioni iniziali
      cinf >> particle >> x >> y >> z; 
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){//se è una simulazione di ISING carico anche gli spins
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk){ // setta a zero gli accumulatori di blocco utilizzati nelle medie di blocco
  ofstream coutf;
  if(blk>0){//scrivo il blocco della simulazione che è stato completato nel file di output
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();//pongo a zero gli accumulatori di blocco
  return;
}

void System :: measure(){ //misuro le proprietà
  _measurement.zeros();//pono a zero il vettore che accumula i valori delle proprietà step per step nel singolo blocco
  vec distance;
  distance.resize(_ndim);
  //accumulatori temporanei delle proprietà da misurare
  double penergy_temp=0.0, dr; 
  double kenergy_temp=0.0; 
  double tenergy_temp=0.0;
  double magnetization=0.0;//non utilizzato in questo esercizio
  double virial=0.0;//non utilizzato in questo esercizio
  if (_measure_penergy or _measure_pressure or _measure_gofr) {//se una di queste 3 proprietà è stata letta dal file properties allora la misuro
    //calcolo per ogni particella l'energia potenziale dovuta all'interazione con le altre particelle
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );//calcolo della distanza della particella i dalla particella j
        if(dr < _r_cut){
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // energia potenziale in modello LJ
          if(_measure_pressure) virial       += 1.0/pow(dr,12) - 0.5/pow(dr,6); // pressione
        }
      }
    }
  }
  // misura della distribuzione delle velocità
  if(_measure_pofv){
    double vmodule=0.0;
    //per ogni particella calcolo il modulo della velocità e riempio il bin in cui tale velocità cade
    for(int i=0;i<_npart;i++){
    vmodule= sqrt(dot(_particle(i).getvelocity(),_particle(i).getvelocity()));
    for(int j=0;j<_n_bins_v;j++){
      if(vmodule>=(j*_bin_size_v)&&vmodule<((j+1)*_bin_size_v)){
       _measurement(_index_pofv+(j))+=1;
       break;
      }
    }
  }
  }
  //energia potenziale finale
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);//vtail è posto a 0 in questo esercizio
    _measurement(_index_penergy) = penergy_temp;
  }
  // energia cinetica
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() );//calcolo dell'energia cinetica media per particella 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // energia totale (cinetica più potenziale)
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;//formula per simulazione MD
    else { //non utilizzata in questo esercizio
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // temperatura
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;//formula delle temperatura in termodinamica
  // pressione
  if (_measure_pressure) _measurement[_index_pressure] = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume;//formula della pressione in termodinamica

  _block_av += _measurement; //aggiorno gli accumulatori di blocco

  return;
}
//in questa funzione calcolo medie di blocco, medie cumulative e incertezza con data blocking di tutte le proprietà da misurare
void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);//medie di blocco
  _global_av  += _average;//somma medie di blocco
  _global_av2 += _average % _average; // somma medie di blocco al quadrato, serve per data blocking
  
  //scrivo su file il blocco, le medie di blocco, le medie cumulative e l'incertezza
  // energia potenziale
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // energia cinetica
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // energia totale
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // temperatura
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //pressione
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //distribuzione delle velocità
  if (_measure_pofv){
    coutf.open("../OUTPUT/pofv.dat",ios::app);
    coutf << "blocco numero " << blk << endl;
    for(int i=0; i<_n_bins_v;i++){
    average  = (_average(_index_pofv+(i)))/(_npart*_bin_size_v);
    sum_average = (_global_av(_index_pofv+(i)))/(_npart*_bin_size_v);
    sum_ave2 = _global_av2(_index_pofv+(i))/(pow(_npart,2)*pow(_bin_size_v,2));
    coutf << setw(12) << (_bin_size_v/2)+i*_bin_size_v
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        }
    coutf.close();
  }
  //percentuale di mosse accettate
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}
//calcolo dell'errore con data blocking
double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}
//restituisce il numero di blocchi della simulazione
int System :: get_nbl(){
  return _nblocks;
}
//restituisce il numero di steps della simulazione
int System :: get_nsteps(){
  return _nsteps;
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
