#include <iostream>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include "random.h"
#include "chromosomes.h"
#include "individuo.h"
#include <fstream>
#include <string>
#include <iomanip>
#include "mpi.h"
#include <sstream> 
#include <vector>

//In questo esercizio risolvo il TSP dei 110 capoluoghi di provincia italiani con 11 core indipendenti che ogni 1000 generazioni si scambiano individui migliori


using namespace std;
using namespace arma;

int main (int argc, char *argv[]){
    //size è il numero totale di processi MPI, di core. Rank la loro etichetta
    int size, rank;
    MPI_Init(&argc,&argv);//inizializza MPI con argc e argv passati dalla linea di comando
    MPI_Comm_size(MPI_COMM_WORLD, &size);//Dice a ciascun processo quanti processi totali partecipano alla comunicazione (cioè quanti core sono attivi).
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//Dice a ciascun processo il proprio ID (rank)
    //inizializzo il generatore di numeri casuali 1 sola volta, quindi al rank = 0
    //nota che lo inizializzo nel rank = 0 ma è inizializzato per tutti i rank che eseguono il codice
    Random rnd;
    if (rank == 0) {
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
                rnd.SetRandom(seed,p1,p2);
    } else cerr << "PROBLEM: Unable to open seed.out" << endl;
    input.close();
    //nel rank 0 apro i file che contengono il nome delle città e le loro longitudini e latitudini
    ifstream city_file("INPUT/prov_ita.txt");
    ifstream coord_file("INPUT/cap_prov_ita.dat");
    //preparo il file dove scriverò:
    //nome città | etichetta numerica | longitudine | latitudine
    ofstream output_file("INPUT/output.dat");

    if (!city_file.is_open() || !coord_file.is_open() || !output_file.is_open() ) {
        cerr << "Errore nell'apertura dei file." << endl;
        return 1;
    }

    string city;
    double lat, lon;
    int index = 1;
    output_file << "#Città" << setw(30) << "Numero" << setw(22) << "Longitudine" << setw(21) << "Latitudine" << endl;

    while (getline(city_file, city) && coord_file >> lon >> lat) {
         output_file << left << setw(22) << city
                << right << setw(12) << index
                << setw(22) << fixed << lon
                << setw(22) << fixed << lat
                << endl;
        index++;
    }
    

    city_file.close();
    coord_file.close();
    output_file.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);// tutti i rank aspettano qua, incluso lo 0
    vector<int> dep(size); // vettore di interi di lunghezza = numero di processi
    iota(dep.begin(), dep.end(), 0); // riempie il vettore con {0, 1, 2, ..., size-1}. il vector rank serve per stabilire chi comunica con chi durante le migrazioni
    //sono etichette dei messaggi inviati e ricevuti. 
    MPI_Status stat;
    MPI_Status stat1;
    //inizializzo chromosomes per ogni rank
    Chromosomes es;
    es.set_rank(rank);
    es.inizialize_rnd();
    es.inizialize_population();
    int Nmigr=1000;//ogni quante generazioni avviene la migrazione
    //genero nuove generazioni
    for(int i=0;i<50000;i++){
        es.new_population();
        if((i+1)%Nmigr==0){//migrazione
            if(rank == 0) {// Solo il rank 0 mischia i rank per decidere le migrazioni 
                rnd.shuffle(dep);
            }
            // Broadcast del vettore dep da rank 0 a tutti gli altri
            MPI_Bcast(dep.data(), size, MPI_INT, 0, MPI_COMM_WORLD);//dep.data() è un puntatore al primo elemento di dep, poi Bcast passerà i successivi size elementi essendo contigui nella memoria
            for(int h=0;h<size;h++){//ogni rank comunica col rank successivo secondo l'ordine di dep
                int sender = dep[h];
                int receiver = dep[(h+1) % size]; //in questo modo va anche da dep[size-1] a dep[0]
                for(int j=0;j<7;j++){//scambio dei primi 7 individui
                    //invio
                    if(rank==sender){
                    mat invio = es.get_individuo(j);//individuo da inviare
                    int dims[2] = {static_cast<int>(invio.n_rows), static_cast<int>(invio.n_cols)};//matrice ha dimensioni rowsXcols
                    MPI_Send(dims, 2, MPI_INT, receiver, 0, MPI_COMM_WORLD);//invio le dimensioni perchè il ricevente non sa a priori le dimensioni della matrice e deve allocare la memoria
                    MPI_Send(invio.memptr(), dims[0]*dims[1], MPI_DOUBLE, receiver, 1, MPI_COMM_WORLD);//invio.memptr() restituisce un puntatore alla memoria contigua dove sono memorizzati i dati della matrice invio
                    }
                    //ricezione
                    if(rank==receiver){
                    int dims[2];
                    MPI_Recv(dims, 2, MPI_INT, sender, 0, MPI_COMM_WORLD, &stat);//riceve le dimensioni della matrice
                    mat ricevo(dims[0], dims[1]);//inizializza la matrice
                    MPI_Recv(ricevo.memptr(), dims[0]*dims[1], MPI_DOUBLE, sender, 1, MPI_COMM_WORLD, &stat1);//riempie la matrice con i dati inviati dal rank sender
                    es.set_individuo(es.get_population()-1-j,ricevo);//ricezione
                    }
                }
                
                }
            }
        
         }
    //ordino la popolazione e salvo il miglior percorso per ogni rank
    es.ascending_order_population();
    ostringstream filename;
    filename << "OUTPUT/best_path_cartesian_coordinates_" << rank << ".dat";
    ofstream coutf;
    coutf.open(filename.str());
    coutf << es.get_individuo(0) << endl;
    coutf.close();
    //ogni processo calcola la sua lunghezza di percorso migliore
    double best_path_lenght=es.distance(0);
    struct {//uso una struct perchè MPI_MINLOC vuole una coppia di dati che non siano separati, cioè variabili diverse. La coppia di dati è il percorso da minimizzare e il rank a cui appartiene
        double value;
        int rank;
    } local, global;//local e global sono due variabii del tipo struct che ho creato, così MPI_MINLOC può lavorare con ndue variabili dello stesso tipo
    //local contiene il dato locale di ogni processo, global conterrà la riduzione globale

    local.value = best_path_lenght;
    local.rank = rank;

    // Riduzione per trovare il minimo valore e chi lo ha trovato
    MPI_Reduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

    //solo il rank = 0 stampa su terminale
    if (rank == 0) {
        cout << "Best path found by rank " << global.rank << " with length " << global.value << endl;
    }

    es.finalize();//finalizzo chromosomes
    MPI_Finalize();//finalizzo MPI
   
    return 0;
}

//nota che la comunicazione è bloccante e va bene così: tutti i processi all'inizio vengono sincronizzati con barrier e 
//la comunicazione avviene a coppie fisse di sender e receiver, uno alla volta. Non c’è rischio di deadlock (blocco reciproco) perché i cicli
//sono organizzati in modo deterministico (cioè per ogni coppia mittente-destinatario è chiaro chi manda e chi riceve).

//MPI_COMM_WORLD in quelle funzioni in cui è chiamato indica che tali funzioni sono chiamate da tutti i rank

//il programma si esegue da terminale con mpirun -np 11 ./main, 11 sono i processi, i rank avviati.
//oppure con mpirun -np 11 --oversubscribe ./main se la macchina di calcolo non ha 11 core indipendenti