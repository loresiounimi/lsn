#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){
   //creo il file su cui caricare i dati
   ofstream outputfile1("randomcontinuumwalk.txt");
   
   if(!outputfile1){
      cout << "Error: Unable to open file randomcontinuumwalk.txt" << endl;
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
  vec posizione(3);//vettore che memorizza e aggiorna la posizione a ogni step MC
  vec r2(N);//vettore che memorizza la media delle medie di blocco
  vec sum(N);//vettore che memorizza il modulo quadro di r
  vec errore(N);//vettore che memorizza l'errore MC
  vec r2_2(N);//vettore che memorizza la media delle medie di blocco al quadrato
   //inizializzo a zero i vettori che fungono da "contenitori"
  r2_2.zeros();
  r2.zeros();
  errore.zeros();

  //angoli che mi servono per muovere il randomwalk con un passo di lunghezza 1 lungo una direzione casuale
  double theta=0.0;
  double phi=0.0;
  for(int y=0;y<N;y++){
      sum.zeros();//pulisco il contenitore sum per avere medie di blocco indipendenti
   for(int j=0;j<l;j++){
      posizione.zeros(); //faccio ripartire il cammino da 0
      for(int i=0;i<n;i++){//simulazione di un cammino di 100 passi, a ogni passo aggiorno la posizione casualmente lungo una direzione qualsiasi
         phi=rnd.Rannyu(0.0,2*M_PI);
         theta=rnd.Rannyu(0.0,M_PI); 
         posizione[0]+=sin(theta)*cos(phi); //x in coord sferiche
         posizione[1]+=sin(theta)*sin(phi); //y in coord sferiche
         posizione[2]+=cos(theta); //z in coord sferiche
         sum[i] += posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2];//calcolo modulo quadro di r
      }
   }
   sum/=l; //calcolo medie di blocco
   r2+=sum; //memorizzo le medie di blocco
   r2_2+=sum%sum; //moltiplicazione tra vec, memorizzo le medie di blocco al quadrato, mi serve per la stima dell'errore MC
}
r2/=n; //calcolo media delle medie di blocco
r2_2/=n;
errore = sqrt((r2_2 - pow(r2, 2)) / N);//calcolo l'errore MC sull'ultimo blocco per ognuno dei 100 passi del random walk, è questo ciò che plotto nel grafico
//carico i dati sul file
for(int i=0;i<n;i++){
   outputfile1 << i+1 << " " << sqrt(r2[i]) << " " << errore[i] << endl;
// Output su file: passo del cammino | media delle medie di blocco | errore statistico sull'ultimo blocco calcolato con data blocking
}
rnd.SaveSeed();//salvo il seed
  cout << 2 << endl;
  return 0;
}

  /*int m=1000000;
  int n=100;
  int l=m/n;
  int step=100;
  double phi=0;
  vector <double> posizione(3);
  //vector <vector<double>> r(2);
  vector <double> sum;
  vector <double> sum_average;
  vector <double> error;
  double dep=0;
  double theta=0;
  //double r;
  for(int y=0;y<n;y++){
   //dep=0;
   for(int j=0;j<l;j++){
      
      for(int i=0;i<step;i++){
         phi=rnd.Rannyu(0.0,2*M_PI);
         theta=rnd.Rannyu(0.0,M_PI); 
         posizione[0]+=sin(theta)*cos(phi); //x in coord sferiche
         posizione[1]+=sin(theta)*sin(phi); //y in coord sferiche
         posizione[2]+=cos(theta); //z in coord sferiche

      }
      //cout << posizione[0] << " " << posizione[1] << " " << posizione[2] << endl;
      dep+=accumulate(posizione.begin(), posizione.end(), 0.0, 
      [](double acc, double val) {
          return acc + val * val;
      }
      );
      posizione.assign(posizione.size(), 0.0);
   }
   sum.push_back(sqrt(dep/static_cast<double>(l)));
   sum_average.push_back(media(sum,y+1));
   error.push_back(errore(sum,y+1));
   //r[0].push_back(sqrt(dep/static_cast<double>(l*(y+1))));
   //r[1].push_back(errore(r[0],y+1));
   //cout << errore(r[0],y+1) << endl;
   outputfile1 << l*(y+1) << " " << sum_average[y] << " " << error[y] << endl;
   //cout << dep/(l*(y+1)) << endl;
}
  //cout << r << " " << posizione[0] << " " << posizione[1] << " " << posizione[2] << endl;
  return 0;
}*/