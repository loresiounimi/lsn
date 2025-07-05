#ifndef __Individuo__
#define __Individuo__
//Individuo è una classe che descrive un individuo della popolazione del TSP. Ogni individuo è un percorso del TSP

#include <iostream>
#include <armadillo>
#include "random.h"
#include <vector>
using namespace std;
using namespace arma;

class Individuo{
private:
int _cities=34; //città da visitare
Random* _rnd; //puntatore a un oggetto Random che in questo caso è il rnd della classe Chromosomes
mat _individuo; //matrice che descrive un individuo

public:
Individuo(); //costruttore
Individuo(const Individuo& other);  // costruttore di copia 
void set_individuo(const mat& m); //setta la matrice che descrive un individuo
void print_individuo(); //stampa l'individuo
mat genera_individuo_iniziale(); //genera l'individuo iniziale per il problema TSP
double costo(); //calcola la funzione costo di un singolo individuo
mat genera_individuo(); //genera un nuovo individuo randomicamente

//mutazioni
void mutazione_semplice(); 
void mutazione_shift(); 
void mutazione_scambio_blocchi();
void mutazione_inversione();

mat rows(int a, int b); //restituisce le righe della matriche individuo dalla a alla b
mat row(int a); //restituisce la riga a della matrice individuo
int get_city_number(int i); //restituisce il numero di città da visitare
void check(); //controlla che il primo individuo abbia la prima città fissa
mat get_individuo(); //retsituisce la matrice delle posizioni delle città che costituiscono un individuo
void set_random(Random& rnd);//imposta l'ogetto Random a cui _rnd punta
void set_cities(int cities);//imposta il numero di città
void carica_individuo_iniziale();//carica da un file l'individuo iniziale
int cities();//restituisce il numero di città
};

#endif // __Individuo__