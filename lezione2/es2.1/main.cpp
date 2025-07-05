#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory> //libreria dei puntatori intelligenti
#include "random.h"
#include "integral.h"
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
   //creo due file su cui caricare i dati della simulazione
   string MC[2]= {"intmedia.txt", "intsampling.txt",};
   ofstream MCwrite[2];
   for(int i=0;i<2;i++){
      MCwrite[i].open(MC[i]);
      if(!MCwrite[i]){
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
  int nmax= 10000; //numero totale di MC steps
  int N=100; //blocchi della simulazione MC
  funzione1 unif; //creo la funzione per calcolare l'integrale MC con campionamento uniforme
  funzione2 sampling; //creo la funzione per calcolare l'integrale MC con campionamento uniforme
  int l=nmax/N; //steps MC per blocco
  vector<unique_ptr<funzione>> a; //unique_ptr<funzione>, è puntatore intelligente a un oggetto di tipo funzione (o a una sua classe derivata). Questo vettore può contenere puntatori a oggetti di tipo funzione, funzione1, funzione2, ecc., purché siano classi derivate da funzione.
    a.push_back(make_unique<funzione1>(unif));  
    a.push_back(make_unique<funzione2>(sampling));
  vector <vector<double>> v(3); //vector che utilizzo per memorizzare x_i, l'integrale MC e l'errore MC
  // v[0] = campioni x_i
  // v[1] = medie degli integrali MC per ogni blocco
  // v[2] = errori statistici calcolati per ogni blocco

  for(int y=0;y<2;y++){ // y = 0: campionamento uniforme, y = 1: importance sampling
   for (int j=0;j<N;j++){
  for (int i=0;i<l;i++){
   if(y==0){
     v[0].push_back(rnd.Rannyu());//memorizzo x_i da campionamento uniforme
   }else{
      v[0].push_back(rnd.funzionesimil());//memorizzo x_i da importance sampling
   } 
   } 
   v[1].push_back(integraleMC(*a[y],v[0],l));//memorizzo l'integrale MC di ogni blocco
   v[2].push_back(errore(v[1], (j+1)));//memorizzo l'errore MC
   
   //scrivo su file i dati ottenuti dalle simulazioni
   MCwrite[y] << sommaElementi(v[1])/(j+1) << " " << l*(j+1) << " " << v[2][j] << endl;
   // Output su file: media cumulativa | passi Monte Carlo totali | errore statistico

   v[0].clear();//pulisco il contenitore delle x_i per avere medie di blocco indipendneti
   }
   //pulisco i contenitori degli integrali MC e dell'errore MC alla fine del ciclo sulle y perchè vengono usati per entrambi i campionamenti 
   v[1].clear();
   v[2].clear();
}

rnd.SaveSeed(); //salvo il seed

  return 0;
}