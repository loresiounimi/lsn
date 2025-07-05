#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"
#include "statistics.h"
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
   //creo due file su cui caricare i dati della simulazione
   string option[2]= {"callprice.txt", "putprice.txt",};
   ofstream optionwrite[2];
   for(int i=0;i<2;i++){
      optionwrite[i].open(option[i]);
      if(!optionwrite[i]){
        cout << "Error: Unable to open the file" << endl;
        return 1;
        } 
     } 
     //inizializzo il generatore di numeri casuali
  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.out");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  int n=100; //numero di blocchi della simulazione Monte-Carlo
  int l=10000; //numero di steps per blocco

  vector<double> dep(4, 0.0);// vettore temporaneo per accumulare i payoff simulati in un blocco (direct/discrete × call/put)

  //inizializzo i parametri della simulazione

  double S0=100.0; //valore iniziale dell'asset
  double K=100.0; //strike
  double T=1.0; //tempo di scadenza
  double r=0.1; //tasso di interesse privo di rischio
  double sigma=0.25; //volatilità

  double dep1;
  //inizializzo le due funzioni di wiener che riproducono il GBM
  wiener S(sigma, r, S0, rnd, 0.0);
  wiener Sdiscrete(sigma,r, S0, rnd, 0.0);

  vector <vector<vector<double>>> price(2); // 1 vector che contiene 2 vector (uno per call, l'altro per put) di vector double ciascuno di dimensione 2 (uno per la media, l'altro per l'incertezza )
  vector <vector<vector<double>>> discreteprice(2); // potevo implementare 1 solo di questi vector, ne implemento 2 per chiarezza del codice

  //inizializzo la dimensione del primo vector contenuto a 2, una per call e una per put
  price[0].resize(2);
  price[1].resize(2);
  discreteprice[0].resize(2);
  discreteprice[1].resize(2);

  //inizia la simulazione Monte-Carlo
  for(int j=0;j<n;j++){
   fill(dep.begin(), dep.end(), 0.0); //setto a 0 i contenitori dep all'inizio di ogni blocco
  for(int i=0;i<l;i++){
   Sdiscrete.setS0(S0); //ripristino a S0 il valore iniziale dell'asset nel campionamento discretizzato perchè il ciclo su t che va da 0 a T utilizza il valore di S(t_i) per calcolare S(t_i+1)
   for(int y=0;y<100;y++){
      Sdiscrete.settprec(static_cast<double>(y)/(100.0));//setto t_i per il calcolo di S(t_{i+1})
      dep1=Sdiscrete.Eval(static_cast<double>(y+1)/(100.0));//calcolo S(t_{i+1})
      Sdiscrete.setS0(dep1);//pongo S0 pari a S(t_i+1) per poter calcolare nel ciclo successivo S(t_{i+2})
     
   }
   //calcolo i payoff attualizzati delle due ozpioni
   dep[0]+=(exp(-r*T)*max(0.0,S.Eval(T)-K));//opzione call, campionamento diretto
   dep[1]+=(exp(-r*T)*max(0.0,K-S.Eval(T)));//opzione put, campionamento diretto
   dep[2]+=(exp(-r*T)*max(0.0,dep1-K));//opzione call, campionamento discretizzato
   dep[3]+=(exp(-r*T)*max(0.0,K-dep1));//opzione put, campionamento discretizzato
}
// calcolo medie di blocco e incertezze (errore statistico) e aggiorno i vettori
price[0][0].push_back(dep[0]/static_cast<double>(l));
price[0][1].push_back(errore(price[0][0],j+1));
price[1][0].push_back(dep[1]/static_cast<double>(l));
price[1][1].push_back(errore(price[1][0],j+1));
discreteprice[0][0].push_back(dep[2]/static_cast<double>(l));
discreteprice[0][1].push_back(errore(discreteprice[0][0],j+1));
discreteprice[1][0].push_back(dep[3]/static_cast<double>(l));
discreteprice[1][1].push_back(errore(discreteprice[1][0],j+1));

//scrivo sui file i dati ottenuti con le simulazioni
optionwrite[0] << l*(j+1) << " " << sommaElementi(price[0][0])/(j+1) << " " << price[0][1][j] << " " << sommaElementi(discreteprice[0][0])/(j+1) << " " << discreteprice[0][1][j] << endl;
optionwrite[1] << l*(j+1) << " " << sommaElementi(price[1][0])/(j+1) << " " << price[1][1][j] << " " << sommaElementi(discreteprice[1][0])/(j+1) << " " << discreteprice[1][1][j] << endl;
//con sommaElementi sommo le componenti del vector e poi le divido per il numero di blocchi a cui è arrivata la simulazione, così ottengo le medie progressive
}

//salvo il seed
 rnd.SaveSeed();
  return 0;
}