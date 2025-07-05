/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//system è una classe in grado di eseguire vari tipi di simulazione molecolare e misurare alcune proprietà medie della simulazione nella sua evoluzione temporale
//in questo esercizio è effettuatta una simulazione del modello di Ising con algoritmi di metropolis e gibbs
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

void System :: move(int i){ // propone una mossa MC alla particella i
  if(_sim_type == 3){ //modello di Ising algoritmo di Gibbs
    //calcolo della probabilità condizionata di flippare lo spin i conoscendo lo spin di i-1 e i+1
    double p_conditional=1/(1+exp(2*_particle(i).getspin()*_beta*_J*(_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() +_H )));
    //se la p è maggiore di un numero estratto casualmente tra 0 e 1 allora flippa lo spin
    if(_rnd.Rannyu() < p_conditional){
      _particle(i).flip();
    } 
    _naccepted++;//aggiorna le mosse accettate
  } else { //algoritmo di metropolis non per il modello di ISING, non usato in questo esercizio
    if(_sim_type == 1){       
      vec shift(_ndim);       
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; 
      }
      _particle(i).translate(shift, _side);  
      if(this->metro(i)){ 
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback();
    } else {                  // Ising 1d con metropolis
      if(this->metro(i)){     //metro restituisce un booleano. true se la mossa è accettata, false se non è accettata
        _particle(i).flip();  //se è accettata, flippa lo spin della particella i
        _naccepted++; //aggiorna le mosse accettate
      }
    }
  }
  return;
}

bool System :: metro(int i){ // algoritmo di metropolis
  bool decision = false; //booleano che stabilisce se accettare o meno la mossa
  //non usato in questo esercizio
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  //metropolis per Ising
  //calcolo differenza di energia tra stato iniziale e stato finale proposto
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  //calcolo della probabilità di accettazione
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //step di accettazione del metropolis
  return decision;//ritorna la decisione. true mossa accettata, false mossa rifiutata
}
//non usato in questo esercizio
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

int System :: pbc(int i){ // applica le pbc per gli spin
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
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl; //Ising con metropolis
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl; //Ising con Gibbs
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
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // setta randomicamente lo spin delle particelle
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
  this->read_configuration();//leggo configurazione iniziale di spins
  if(_sim_type==0) this->initialize_velocities();//inizializzo le velocità se è MD Lennard-Jones
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
  } else {//se non è presente il file con xold, il codice assegna le velocità casualmente ponendole uguali a +- velocity, in questo modo parto da uno stato iniziale a bassa entropia per quanto riguarda le velocità inziali
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    double direction;
    double verse;
    double velocity=sqrt(3*_temp); //velocità impostata a +- radice di 3*T per rispettare il legame tra energia cinetica media del gas e temperatura
    sumv.zeros();
    //assegno le velocità
    for (int i=0; i<_npart; i++){
      direction=static_cast<int>(_rnd.Rannyu(1.0,4.0));// 1 direzione x; 2 direzione y; 3 direzione z
      verse=static_cast<int>(_rnd.Rannyu(1.0,3.0)); // 1 velocità positiva; 2 velocità negativa
      if(verse==1){
        if(direction==1){vx(i)=velocity;}  
        if(direction==2){vy(i)=velocity;}
        if(direction==3){vz(i)=velocity;}
     }else if(verse==2){
        if(direction==1){vx(i)=-velocity;}
        if(direction==2){vy(i)=-velocity;}
        if(direction==3){vz(i)=-velocity;}
     }
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
      } else if( property == "MAGNETIZATION" ){//magnetizzazione
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){//calore specifico
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){//suscettività
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
  } else {//per una simulazione ISING
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
  cinf.open("../INPUT/CONFIG/config.ising");//apro il file
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
  } else cerr << "PROBLEM: Unable to open INPUT file config.ising"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){//se è una simulazione di ISING memorizzo gli spins
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
    if(blk==1){coutf << "Temperatura= " << _temp << endl;}//per il primo blocco di ogni simulazione scrivo su output.dat la temperatura a cui la simulazione è eseguita
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();//pongo a zero gli accumulatori di blocco
  return;
}


void System :: measure(){ // misuro le proprietà
  _measurement.zeros();//pono a zero il vettore che accumula i valori delle proprietà step per step nel singolo blocco
  vec distance;
  distance.resize(_ndim);
  //accumulatori temporanei delle proprietà da misurare
  double penergy_temp=0.0, dr; 
  double kenergy_temp=0.0; 
  double tenergy_temp=0.0;
//accumulatore del quadrato dell'energia media per particella, serve per il calcolo del cv
  double tenergy_temp2=0.0;
//accumulatore degli spins, serve per il calcolo della suscettività
  double s_temp=0.0;
  double magnetization=0.0;
  double virial=0.0;//non usato in questo esercizio
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
  // energia totale 
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;//formula per simulazione MD
    else { //energia media per particella in una simulazione ISING in cui c'è interazione tra spins vicini, in questa simulazione 1D l'interazione è solo tra primi vicini
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
  //magnetizzazione
  if (_measure_magnet){
    if (_sim_type > 1)  { 
      for (int i=0; i<_npart; i++){
        magnetization += double(_particle(i).getspin());
      }
      _measurement(_index_magnet) = magnetization/double(_npart) ; 
    }
  }
  //calore specifico
  if (_measure_cv){//in questo caso viene memorizzata solo l'energia media per particella ^2, il calcolo finale di cv è in averages
    if (_sim_type > 1)  { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp2 += pow(- _J * s_i * s_j - 0.5 * _H * (s_i + s_j),2);
      }
      tenergy_temp2 /= double(_npart);
      _measurement(_index_cv) = tenergy_temp2 ; 
    }
  }
  //suscettività
  if (_measure_chi){
    if (_sim_type > 1)  { 
      for (int i=0; i<_npart; i++){
        s_temp += double(_particle(i).getspin());
      }
      _measurement(_index_chi) = _beta*(pow(s_temp,2)/double(_npart)) ; 
    }
  }

  _block_av += _measurement; //aggiorno gli accumulatori di blocco

  return;
}
//in questa funzione calcolo medie di blocco, medie cumulative e incertezza con data blocking di tutte le proprietà da misurare
void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);//medie di blocco
  //calcolo corretto di cv per le medie di blocco
  if(_measure_cv){
  double dep=_average(_index_cv);
  _average(_index_cv)  = pow(_beta,2)*(dep-pow(_average(_index_tenergy),2)); //correct average for specific heat
  }
   _global_av  += _average;//somma medie di blocco
  _global_av2 += _average % _average; // somma medie di blocco al quadrato, serve per data blocking

  //scrivo su file il blocco, le medie di blocco, le medie cumulative e l'incertezza
  //nelle grandezze misurate in questo esercizio di simulazione ISING, al primo blocco scrivo
  //anche la temperatura a cui avviene la simulazione
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
  //energia cinetica
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
  //energia totale
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    if(blk==1){coutf << "# Temperatura= " << _temp << endl;}
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //temperatura
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
  //magnetizzazione
  if (_measure_magnet){
    coutf.open("../OUTPUT/magnetization.dat",ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    if(blk==1){coutf << "# Temperatura= " << _temp << endl;}
    coutf << setw(20) << blk
          << setw(20) << average
          << setw(20) << sum_average/double(blk)
          << setw(20) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //calore specifico
  if (_measure_cv){
    coutf.open("../OUTPUT/specific_heat.dat",ios::app);
    average  = _average(_index_cv);
    sum_average = _global_av(_index_cv);
    sum_ave2 = _global_av2(_index_cv);
    if(blk==1){coutf << "# Temperatura= " << _temp << endl;}
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //suscettività
  if (_measure_chi){
    coutf.open("../OUTPUT/susceptibility.dat",ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    if(blk==1){coutf << "# Temperatura= " << _temp << endl;}
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //percentuale di mosse accettate
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(blk==1){coutf << "Temperatura= " << _temp << endl;}
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
//imposta il valore della temperatura
void System :: set_temperature(double t){
  _temp=t;
}
//imposta il valore di beta
void System :: set_beta(double t){
  _beta=1.0/t;
}
//restituisce la temperatura
double System ::  get_temperature(){
  return _temp;
}
//azzera gli accumulatori di medie cumulative ed errore
void System :: global_reset(){
  _global_av.zeros();
  _global_av2.zeros();
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
