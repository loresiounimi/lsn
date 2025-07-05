#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){
   //creo il file su cui caricare i dati
   ofstream outputfile1("randomdiscretewalk.txt");
   
   if(!outputfile1){
      cout << "Error: Unable to open file randomdiscretewalk.txt" << endl;
      return 1;
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
  int n=100; //numero di passi di ogni simulazione
  int N=100; //numero di blocchi della simulazione
  int l=10000; //simulazioni per blocco
  vec posizione(3); //vettore che memorizza e aggiorna la posizione a ogni step MC
  vec r2(N); //vettore che memorizza la media delle medie di blocco
  vec sum(N); //vettore che memorizza il modulo quadro di r
  vec errore(N); //vettore che memorizza l'errore MC
  vec r2_2(N);//vettore che memorizza la media delle medie di blocco al quadrato
  //inizializzo a zero i vettori che fungono da "contenitori"
  r2_2.zeros();
  r2.zeros();
  errore.zeros();

  int direzione=0;
  int verso=0;
  for(int y=0;y<N;y++){
      sum.zeros(); //pulisco il contenitore sum per avere medie di blocco indipendenti
   for(int j=0;j<l;j++){
      posizione.zeros();//faccio ripartire il cammino da 0
      for(int i=0;i<n;i++){ //simulazione di un cammino di 100 passi, a ogni passo aggiorno la posizione casualmente lungo un reticolo discreto di lato 1
         direzione=static_cast<int>(rnd.Rannyu(1.0,4.0));// 1 direzione x; 2 direzione y; 3 direzione z
         verso=static_cast<int>(rnd.Rannyu(1.0,3.0)); // 1 passo in avanti; 2 passo indietro
         if(verso==1){
            if(direzione==1){posizione[0]+=1.0;}
            if(direzione==2){posizione[1]+=1.0;}
            if(direzione==3){posizione[2]+=1.0;}
         }else if(verso==2){
            if(direzione==1){posizione[0]-=1.0;}
            if(direzione==2){posizione[1]-=1.0;}
            if(direzione==3){posizione[2]-=1.0;}
         }
         sum[i] += posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]; //calcolo modulo quadro di r
      }
   }
   sum/=l; //calcolo medie di blocco
   r2+=sum; //memorizzo le medie di blocco
   r2_2+=sum%sum; //moltiplicazione tra vec, memorizzo le medie di blocco al quadrato, mi serve per la stima dell'errore MC
}
r2/=n; //calcolo media delle medie di blocco
r2_2/=n;
errore = arma :: sqrt((r2_2 - pow(r2, 2)) / N); //calcolo l'errore MC sull'ultimo blocco per ognuno dei 100 passi del random walk, è questo ciò che plotto nel grafico
//carico i dati sul file
for(int i=0;i<n;i++){
   outputfile1 << i+1 << " " << sqrt(r2[i]) << " " << errore[i] << endl;
   // Output su file: passo del cammino | media delle medie di blocco | errore statistico sull'ultimo blocco calcolato con data blocking
}
rnd.SaveSeed();//salvo il seed
  
  return 0;
}