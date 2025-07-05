/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//system è una classe preparata per simulazioni di dinamica molecolare e modificata
//per il calcolo di <H>
#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // esegue un passo della simulazione
  for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // esegue un passo di simulazione MC su una particella scelta a caso
  _nattempts += _npart; //aggiorna il numero di tentativi eseguiti sul sistema
  return;
}

void System :: move(int i){ // propone una mossa MC alla particella i
    //algoritmo di metropolis per simulazione Montecarlo    
      //propone la mossa MC   
      vec shift(_ndim);       
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; 
      }
      _particle(i).translate(shift);  
      //se la mossa viene accettata, la posizione viene aggiornata. Se la mossa viene rifiutata la particella torna alla posizione iniziale
      if(this->metro(i)){ 
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback();
    
  return;
}

bool System :: metro(int i){ // algoritmo di metropolis
  bool decision = false; //booleano che stabilisce se accettare o meno la mossa
  //calcolo probabilità di accettazione per simulazione MC
  double acceptance;
  acceptance = pow(fabs(this->PSI(i,true)),2)/pow(fabs(this->PSI(i,false)),2);
  if(_rnd.Rannyu() < acceptance ) decision = true; //step di accettazione del metropolis
  return decision;//ritorna la decisione. true mossa accettata, false mossa rifiutata
}

//funzione che calcola PSI(x)
double System :: PSI(int i, bool xnew){
  return exp(-pow((_particle(i).getposition(0,xnew)-_media),2)/(2*pow(_sigma,2)))+exp(-pow((_particle(i).getposition(0,xnew)+_media),2)/(2*pow(_sigma,2)));
}

//funzione che calcola il modulo quadro di PSI(x)
double System :: PSI_seconda(int i, bool xnew){
  return ((pow((_particle(i).getposition(0,xnew)-_media),2))/pow(_sigma,4)-1/pow(_sigma,2))*exp(-pow((_particle(i).getposition(0,xnew)-_media),2)/(2*pow(_sigma,2)))+((pow((_particle(i).getposition(0,xnew)+_media),2))/pow(_sigma,4)-1/pow(_sigma,2))*exp(-pow((_particle(i).getposition(0,xnew)+_media),2)/(2*pow(_sigma,2)));
}

void System :: initialize(){ // inizializza il sistema in base al contenuto dei file di input che si trovano nella cartella INPUT

  int p1, p2; // inizializza il rnd
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; 
  ifstream Seed("../INPUT/seed.out");
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
    if( property == "RESTART" ){ //restart dice al sistema se generare valori per le velocità e dunque per le posizioni precdenti oppure se leggere le posizioni precedenti da un file già esistente
      input >> _restart;
    } else if( property == "TEMP" ){ //memorizza la temperatura e il suo inverso beta
      input >> _temp;
      _beta = 1.0/_temp;
    } else if( property == "NPART" ){//memorizza il numero di particelle
      input >> _npart;
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "DELTA" ){//passo della simulazione MC
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){//numero di blocchi della simulazione
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){//numero di steps per blocco
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "MEDIA" ){ //valore di mu
      input >> _media;
      coutf << "MEDIA= " << _media << endl;
    } else if( property == "SIGMA" ){ //valore di sigma
      input >> _sigma;
      coutf << "SIGMA= " << _sigma << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;//scrive sul file di output che la letture dei dati di input è andata abuon fine
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();//leggo la configurazione iniziale
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_properties(){ // inizializza i data membri utilizzati per la misura delle proprietà del sistema
  //ogni proprietà da misurare ha il suo indice, serve per poter accedere ai vettori che svolgono le operazioni di media e data blocking
  //in questo codice, l'unica proprietà da misurare è <H>
  string property;
  int index_property = 0;
  _nprop = 0;//inizialmente il numero di proprietà è pari a zero

  _measure_H = false;
  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
     if( property == "H" ){//valor medio di H
        ofstream coutpr("../OUTPUT/<H>.dat");//apro il file e carico i dati
        coutpr << "#     BLOCK:   ACTUAL_H:     H_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_H = true;
        _index_H = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
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
    coutf.open("../OUTPUT/CONFIG/config.xyz");//apro il file
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){//carico le posizioni finali
              coutf << setprecision(17) << _particle(i).getposition(0,true) << endl; // x

      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");//carico le posizioni precedenti alle posizioni finali
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << setprecision(17) << _particle(i).getposition(0,false) << endl;// x
        }    
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
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
    double x;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;

    for(int i=0; i<_npart; i++){
      cinf >> particle >> x; //carico le posizioni iniziali
      _particle(i).setposition(0, x);
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.ising"<< endl;
  cinf.close();
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

void System :: measure(){ // misura delle proprietà
  _measurement.zeros();//pono a zero il vettore che accumula i valori delle proprietà step per step nel singolo blocco
  //accumulatore temporaneo di H
  double H_temp= 0.0;
  //misura di H
  if (_measure_H){
       for (int i=0; i<_npart; i++){
        //calcolo del valor medio di H
        H_temp += -0.5*(this->PSI_seconda(i,true)/this->PSI(i,true))+pow(_particle(i).getposition(0,true),4)-(5.0/2.0)*pow(_particle(i).getposition(0,true),2);
      }
       _measurement(_index_H) = H_temp/double(_npart);
       }
  _block_av += _measurement; //aggiorno gli accumulatori di blocco

  return;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);//medie di blocco
   _global_av  += _average;//somma medie di blocco
  _global_av2 += _average % _average; // somma medie di blocco al quadrato, serve per data blocking 
  //<H>
  if (_measure_H){
    coutf.open("../OUTPUT/<H>.dat",ios::app);
    average  = _average(_index_H);
    sum_average = _global_av(_index_H);
    sum_ave2 = _global_av2(_index_H);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
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
//imposta il valore della temperatura
void System :: set_temperature(double t){
  _temp=t;
}
//imposta il valore di beta
void System :: set_beta(double t){
  _beta=1.0/t;
}
//imposta il valore di mu
void System :: set_media(double media){
  _media=media;
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "new media = " << _media << endl;
  coutf.close();
}
//imposta il valore di sigma
void System :: set_sigma(double sigma){
  _sigma=sigma;
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "new sigma = " << _sigma << endl;
  coutf.close();
}
//restituisce il valore di mu
double System ::  get_media(){
  return _media;
}
//restituisce il valore di sigma
double System ::  get_sigma(){
  
  return _sigma;
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
//restituisce la posizione 1D (per convenzione x) della particella i 
double System :: get_x(int i){
  return _particle(i).getposition(0,true);
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
