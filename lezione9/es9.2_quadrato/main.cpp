#include <iostream>
#include "random.h"
#include "chromosomes.h"
#include "individuo.h"
#include <fstream>

//in questo esercizio si risolve il TSP di 34 città posizionate entro un quadrato di lato 2

using namespace std;

int main (int argc, char *argv[]){
  //creo un oggetto Chromosomes e inizializzo la popolazione
  Chromosomes es;
  es.inizialize_rnd();
  es.inizialize_population();

  cout << es.distance(0) << endl;
  for(int i=0;i<10000;i++)es.new_population();//genero nuove popolazioni 
  cout << es.distance(0) << endl;
  //salvo su un file le posizioni delle città del miglior individuo
  ofstream coutf;
  coutf.open("OUTPUT/best_path_cartesian_coordinates.dat");
  coutf << es.get_individuo(0) << endl;
  coutf.close();
  es.finalize();//salvo il seed


  
    return 0;
    }