#include "chromosomes.h"
#include "random.h"
#include <iomanip>

using namespace std;
using namespace arma;

//inizializza il generatore di numeri casuali
void Chromosomes :: inizialize_rnd(){
    int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.out");
  if (input.is_open()){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           _rnd.SetRandom(seed,p1,p2);
  } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  input.close();
  }

//genera la prima popolazione
    void Chromosomes :: inizialize_population(){
        _individuo.set_random(_rnd);//imposta il rnd di individuo
        _individuo.genera_individuo_iniziale();//genera l'individuo iniziale
        _popolazione.set_size(_population);//imposta le dimensioni dell'array di individui
        //inizializzo i file su cui scrivo l'andamento della miglior stima di L2(x) e di <L2(x)> calcolato su metà popolazione
        ofstream coutf;
        coutf.open("OUTPUT/best_path.dat");
        ofstream coutf1;
        coutf1.open("OUTPUT/best_path_half_population.dat");
        //inizializzo la prima popolazione
        for(int i=0;i<_population;i++){
            _individuo.set_individuo(_individuo.genera_individuo());
            _popolazione[i]=_individuo;
            }
           _ngeneration++;//aggiorno il numero di generazioni
           _popolazione_vecchia=_popolazione;
           this->ascending_order_population();//ordina la popolazione in modo crescente sulla base della funzione costo dei singoli individui
           //scrivo sul file la funzione costo del milgior individuo della prima generazione
           coutf << "#Generation" << setw(20) << "L1(x) best" << endl;
           coutf << _ngeneration << setw(29) << _popolazione[0].costo() << endl;
           //calcolo <L2(x)> della prima generazione
           double dep=0.0;
           double r = static_cast<double>(_population)/2.0;
           for(int i=0;i<static_cast<int>(r);i++){
              dep+=_popolazione[i].costo();
            }
          //scrivo sul file
          coutf1 << "#Generation" << setw(20) << "<L1(x)> best" << endl;
          coutf1 << _ngeneration << setw(29) << dep/r << endl;
          coutf.close();
          coutf1.close();
        }

Chromosomes::Chromosomes() {
    //costruttore vuoto
}
//restituisce la funzione costo di un individuo della popolazione
double Chromosomes :: distance(int individuo){
    return _popolazione[individuo].costo();
}
//stampa un individuo della popolazione
void Chromosomes :: print_individuo(int individuo){
  _popolazione[individuo].print_individuo();
}
//ordina la popolazione in modo crescente sulla base della funzione costo
void Chromosomes :: ascending_order_population(){
  for(int i=0;i<_population-1;i++){
      for(int j=i+1;j<_population;j++){
        if(_popolazione[j].costo()<_popolazione[i].costo()){
                Individuo temp = _popolazione[i];
                _popolazione[i] = _popolazione[j];
                _popolazione[j] = temp;
          }
        }
    }     
    _popolazione_vecchia=_popolazione;
  }
//imposta il numero di individui della popolazione
void Chromosomes :: set_population(int a){
  _population=a;
  }
//restituisce il numero di individui della popolazione
int Chromosomes :: get_population(){
  return _population;
  }

//operatore di crossover, genero due figli da due genitori
void Chromosomes::crossover(Individuo& padre, Individuo& madre){
    //posizione da cui prendere i blocchi iniziali di padre e madre
    int cut = static_cast<int>(_rnd.Rannyu(1, _cities - 1));

    // Prendo il blocco iniziale da padre e madre
    mat bloccoPadre = padre.rows(0, cut);
    mat bloccoMadre = madre.rows(0, cut);

    vector<int> citta_in_blocco_padre;
    vector<int> citta_in_blocco_madre;

    // Salvo le città già prese nel blocco iniziale (da 0 a cut)
    for(int i = 0; i <= cut; i++){
        citta_in_blocco_padre.push_back(padre.get_city_number(i));
        citta_in_blocco_madre.push_back(madre.get_city_number(i));
    }

    // Completo figlio1: aggiungo le città di madre che non sono già in bloccoPadre seguendo l'ordine in cui compaiono in madre
    for(int i = 0; i < _cities; i++){
        int city = madre.get_city_number(i);
        if(find(citta_in_blocco_padre.begin(), citta_in_blocco_padre.end(), city) == citta_in_blocco_padre.end()){ //find è una funzione che cerca city utilizzando gli iteratori .begin e .end. se non trova city restituircse .end() cioè arriva alla fine del vector considerato
            //e aggiunge la città mancante in blocco padre
            bloccoPadre = join_cols(bloccoPadre, madre.row(i));
            //aggiorna le città in blocco padre
            citta_in_blocco_padre.push_back(city);
        }
    }

    // Completo figlio2: aggiungo le città di padre che non sono già in bloccoMadre seguendo l'ordine in cui compaiono in padre. stessa operazione effettuata per il
    //figlio1 ma con madre e padre invertiti
    for(int i = 0; i < _cities; i++){
        int city = padre.get_city_number(i);
        if(find(citta_in_blocco_madre.begin(), citta_in_blocco_madre.end(), city) == citta_in_blocco_madre.end()){
            bloccoMadre = join_cols(bloccoMadre, padre.row(i));
            citta_in_blocco_madre.push_back(city);
        }
    }

    // Imposto i nuovi individui
    padre.set_individuo(bloccoPadre);
    madre.set_individuo(bloccoMadre);
}
//mischia casualmente gli individui nella popolazione
void Chromosomes::shuffle_population() {
    for (int i = _population - 1; i > 0; --i) {
        int j = static_cast<int>(_rnd.Rannyu(0, i + 1)); // intero tra 0 e i
        swap(_popolazione(i), _popolazione(j));
    }
}
//genera una nuova popolazione attraverso mutazioni e crossover
void Chromosomes :: new_population(){
  this->ascending_order_population();//ordino la popolazione in ordine crescente di funzione costo
  //preparo i file su cui scrivere L2(x) del miglior individuo e <L2(x)>
  ofstream coutf;
  ofstream coutf1;
  //genero una nuova popolazione applicando elitismo, salvando cioè i migliori 5 individui in questo caso.
  //inizialmente, senza mutazioni, la nuova popolazione è generata selezionando individui dalla vecchia popolazione attarverso l'operazione di selezione
  //che favorisce i migliori individui 
  for(int i=_elitism;i<_population;i++){
     _popolazione[i]=_popolazione_vecchia[this->selection_operator()]; 
    }
    //applico il crossover con probabilità del 50%
  for(int i=_elitism;i<_population-2;i+=2){
    if(_rnd.Rannyu(0.0,1.0)>0.5)this->crossover(_popolazione[i],_popolazione[i+1]);
}
//due possibili applicazioni di mutazioni sulla popolazione

//1
//applico mutazioni con probabilità del 10%
/*for(int i=_elitism;i<_population;i++){
  if(_rnd.Rannyu(0.0,1.0)>0.9){ //noto che con maggiore probabilità di mutazione ottengo un risultato migliore di L(x) con un numero minore di generazioni. Ne risente <L1(x)> che ha maggiori fluttuazioni.
    int mutation = static_cast<int>(_rnd.Rannyu(1.0,5.0));
    if(mutation==1)_popolazione[i].mutazione_semplice();
    if(mutation==2)_popolazione[i].mutazione_shift();
    if(mutation==3)_popolazione[i].mutazione_scambio_blocchi();
    if(mutation==4)_popolazione[i].mutazione_inversione();
    }
  }*/

  //2
  //applico mutazioni con probabilità del 20% ciascuna
  for(int i=_elitism;i<_population;i++){
  if(_rnd.Rannyu(0.0,1.0)>0.8){ 
    _popolazione[i].mutazione_semplice();
    }
  }
  for(int i=_elitism;i<_population;i++){
    if(_rnd.Rannyu(0.0,1.0)>0.8){
    _popolazione[i].mutazione_shift();
    }
    }
  
  for(int i=_elitism;i<_population;i++){
  if(_rnd.Rannyu(0.0,1.0)>0.8){ 
    _popolazione[i].mutazione_scambio_blocchi();
    }
  }
  for(int i=_elitism;i<_population;i++){
  if(_rnd.Rannyu(0.0,1.0)>0.8){
    _popolazione[i].mutazione_inversione();
    }
  }
  //aggiorno il numero di generazioni
  _ngeneration++;
  _popolazione_vecchia=_popolazione;
  this->ascending_order_population();//ordino la popolazione in ordine crescente della funzione costo
  //apro i file su cui scrivere L2(x) del miglior individuo e <L2(x)> e carico i dati
  coutf.open("OUTPUT/best_path.dat",ios::app);
  coutf << _ngeneration << setw(29) << _popolazione[0].costo() << endl;
  coutf1.open("OUTPUT/best_path_half_population.dat",ios::app);
  double dep=0.0;
  double r=static_cast<double>(_population)/2.0;
  for(int i=0;i<static_cast<int>(r);i++){
      dep+=_popolazione[i].costo();
    }
    coutf1 << _ngeneration << setw(29) << dep/r << endl;
    coutf.close();
    coutf1.close();
  

}
//applica una delle 4 mutazioni di individuo sull'individuo i
void Chromosomes :: mutation (int i){
    int mutation = static_cast<int>(_rnd.Rannyu(1.0,5.0));
    if(mutation==1)_popolazione[i].mutazione_semplice();
    if(mutation==2)_popolazione[i].mutazione_shift();
    if(mutation==3)_popolazione[i].mutazione_scambio_blocchi();
    if(mutation==4)_popolazione[i].mutazione_inversione();
  }
//operatore di selezione che favorisce la selezione di individui con funzione costo più bassa in una popolazione ordinata in modo crescente secondo la funzioen costo stessa
int Chromosomes :: selection_operator(){
    double r = _rnd.Rannyu(0.0, 1.0);
    int j = static_cast<int>(static_cast<double>(_population) * pow(r, 1.7));
    return j;
  }

    //restituisce un individuo della popolazione
mat Chromosomes :: get_individuo(int Individuo){
  return _popolazione(Individuo).get_individuo();
  }
//imposta l'elitismo
void Chromosomes :: set_elitism(int a){
  _elitism=a;
  }
//restituisce l'elitismo
int Chromosomes :: get_elitism(){
  return _elitism;
  }
//finalizza la popolazione salvando il seed
void Chromosomes :: finalize(){
  _rnd.SaveSeed();
  }