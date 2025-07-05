#ifndef __Chromosomes__
#define __Chromosomes__
//questa classe contiene una popolazione di oggetti individuo e genera nuove popolazioni secondo l'algoritmo genetico per risolvere il TSP

#include <iostream>
#include <armadillo>
#include "random.h"
#include "individuo.h"
#include <vector>

using namespace std;
using namespace arma;

class Chromosomes{
private:
int _cities=11; //città di un individuo
int _population=10; //popolazione
Random _rnd; //generatore di numeri casuali
Individuo _individuo; //oggetto individuo
field <Individuo> _popolazione; //contenitore di individui. è una popolazione
field <Individuo> _popolazione_vecchia;//popolazione prima delle mutazioni 
int _elitism=5; //elitismo
int _ngeneration=0; //numero di generazioni
int _rank=0; //etichetta del core

public:
Chromosomes();//costruttore vuoto
void set_rank(int rank);//imposta l'etichetta del core
void inizialize_rnd();//inizializza il rnd
void set_cities(int cities);//imposta il numero di città
double distance(int individuo);//calcola la funzione costo di un individuo
void inizialize_population();//inizializza la popolazione iniziale
void print_individuo(int individuo);//stampa un individuo della popolazione
void ascending_order_population();//ordina la popolazione in ordine crescente sulla base della funzione costo
void set_population(int a);//imposta il numero di individui nella popolazione
int get_population();//restituisce il numero di inividui nella popolazione
void new_population();//genera una nuova popolazione
void crossover(Individuo& padre, Individuo& madre);//operatore di crossover
void shuffle_population();//mischia la popolazione
void mutation(int i);//genera una mutazione
int selection_operator();//operatore di selezione sugli individui
mat get_individuo(int Individuo);//restituisce la matrice di posizioni di un individuo
void set_elitism(int a);//imposta l'elitismo
int get_elitism();//restituisce l'elitismo
void finalize();//finalizza la popolazione
void set_individuo(int i, mat a);//imposta l'individuo i della popolazione
};

#endif // __Chromosomes__