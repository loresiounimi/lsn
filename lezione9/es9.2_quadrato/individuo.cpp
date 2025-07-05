#include "individuo.h"
#include "random.h"

using namespace std;
using namespace arma;

//matrice che genera l'individuo iniziale secondo le condizioni del problema TSP
mat Individuo :: genera_individuo_iniziale(){
    double x=0.0;
    double y=0.0;

    _individuo.set_size(_cities, 3); 
    //numero etichetta della città  || posizione x  || posizione y

    //generazione primo individuo con le città poste entro un quadrato di lato 1
    for(int i=0;i<_cities;i++){
        x=_rnd->Rannyu(-1.0,1.0);
        y=_rnd->Rannyu(-1.0,1.0);
        _individuo(i,0)=i+1;
        _individuo(i,1)=x;
        _individuo(i,2)=y;
        }
    this->check();// controllo sulla prima città
    return _individuo;
    }

    //a partire dal primo inividuo, genero un altro individuo mischiando randomicamente l'ordine delle città
mat Individuo :: genera_individuo(){
    
    mat nuovo_individuo;
    nuovo_individuo=_individuo;
    uvec shuffled = shuffle(regspace<uvec>(1, _cities - 1));//uvec è un vettore di interi. shuffle mescola il vettore di interi creato con indici da 1 a 33. l'indice 0 è fisso essendo la posizione iniziale
    //sostituisco le righe da 1 a 33 (1 è la seconda riga) di nuovo individuo con le righe di individuo in ordine shuffled
    nuovo_individuo.rows(1, _cities - 1) = _individuo.rows(shuffled);
    if(nuovo_individuo(0,0)!=1){
        cout << "errore nella generazione dell'individuo!" << endl;
        return nuovo_individuo.zeros();
    }
    //cambio la matrice individuo
    _individuo=nuovo_individuo;
    this->check();//controllo sulla prima città
    return _individuo;
    }
    //costruttore di copia di Individuo
Individuo :: Individuo(const Individuo& other) {
    _rnd=other._rnd;    // Copia il puntatore, cioè fa sì che punti allo stesso Random dell'individuo copiato
    _cities = other._cities;  // Copia il numero di città
    _individuo = other._individuo;    // Copia la matrice 
    this->check(); //controllo sulla prima città
}
    
//imposto la matrice delle posizioni di un individuo
void Individuo :: set_individuo(const mat& m) {
    _individuo = m;
    _cities = 34; 
    this->check();

}

//stampa della matrice individuo
void Individuo :: print_individuo(){
    _individuo.print();
    }
//costruttore vuoto
Individuo :: Individuo() {
    _cities = 34;
    _individuo.set_size(_cities, 3);
}
//funzione costo di individuo basata su L1(x)
double Individuo :: costo() {
    double costo=0.0;
    for(int i=0;i<_cities-1;i++){
        costo+=sqrt(pow(_individuo(i,1)-_individuo(i+1,1),2)+pow(_individuo(i,2)-_individuo(i+1,2),2));
        }
        costo+=sqrt(pow(_individuo(0,1)-_individuo(_cities-1,1),2)+pow(_individuo(0,2)-_individuo(_cities-1,2),2));  
    return costo;
}
//scambio di due città
void Individuo :: mutazione_semplice(){
    double change1 = static_cast<double>(_rnd->Rannyu(1.0, 34.0));
    int change2 = static_cast<int>(_rnd->Rannyu(1.0, 34.0));
    if (change1 != change2) {
        _individuo.swap_rows(3, change2);// swap rows scambio le due righe selezionate da individuo
    }
    this->check();
}
//shift di n posizioni di m città contigue
void Individuo :: mutazione_shift(){
    while(true){
    int max_m = _cities - 2;  // perché la prima città è fissa
    if (max_m < 2) return; // troppe poche città per fare lo shift

    int m = static_cast<int>(_rnd->Rannyu(1, max_m)); // lunghezza blocco
    int shift = static_cast<int>(_rnd->Rannyu(1, max_m)); // di quante posizioni shiftare
    int start = static_cast<int>(_rnd->Rannyu(1, _cities - m)); // da dove parte il blocco 

    if (start + m + shift >= _cities) continue; // non ci sta il blocco dopo lo shift, riprova con nuovi valori

    // Estraggo il blocco da spostare
    mat blocco = _individuo.rows(start, start + m - 1);
    

    // Cancello il blocco
    mat restante = join_cols(//concatena due matrici verticalmente, cioè una sopra l’altra.
        _individuo.rows(0, start - 1),
        _individuo.rows(start + m, _cities - 1)
    );

    // Ricostruisco l’individuo inserendo il blocco in posizione nuova
    mat nuovo_individuo;

    int insert_pos = start + shift;

    nuovo_individuo = join_cols(
        restante.rows(0, insert_pos - 1),
        blocco
    );
    nuovo_individuo = join_cols(
        nuovo_individuo,
        restante.rows(insert_pos, restante.n_rows - 1)
    );

    _individuo = nuovo_individuo;
    break;
    }
    this->check();
}
//scambio di due blocchi di città in due posizioni diverse
void Individuo :: mutazione_scambio_blocchi(){
    while(true){
        // Calcolo la dimensione massima di un blocco da scambiare (al massimo metà del percorso, meno 2)
        int max_m = (_cities - 2)/2;
        if (max_m < 2) return;
        // Estrazione casuale della lunghezza del blocco da scambiare (tra 1 e max_m)
        int m = static_cast<int>(_rnd->Rannyu(1, max_m));
        // Estrazione casuale della posizione di partenza del primo blocco (tra 1 e _cities - m)
        int start = static_cast<int>(_rnd->Rannyu(1, _cities - m));
        // Estrazione casuale della posizione di partenza del secondo blocco (dopo il primo blocco + m)
        int start2 = static_cast<int>(_rnd->Rannyu(start + m + 1, _cities));
       // Controllo che il secondo blocco non superi il numero totale di città
        if (start2 + m >= _cities) continue;  // Se sì, ricomincio il ciclo
        
        // Estragg1 i blocchi per riordinarli
        mat blocco = _individuo.rows(start, start + m - 1);
        mat blocco1 = _individuo.rows(start2, start2 + m - 1);
        mat blocco2 = _individuo.rows(0, start - 1);
        mat blocco3 = _individuo.rows(start2+m, _cities - 1);
        mat blocco5 = _individuo.rows(start+m, start2 - 1);

        // Ricostruisco la matrice 'nuovo_individuo' concatenando i blocchi:
        // prima blocco2 (parte iniziale), poi blocco1 (secondo blocco spostato in posizione del primo),
        // poi blocco5 (blocco intermedio rimasto invariato), poi blocco (primo blocco spostato in posizione del secondo),
        // infine blocco3 (parte finale)
        mat nuovo_individuo;

        nuovo_individuo = join_cols(join_cols(join_cols(join_cols(blocco2,blocco1),blocco5),blocco),blocco3);

        _individuo = nuovo_individuo;
        break;
        }
        this->check();//controllo sulla prima città
    }
    //inversione dell'ordine delle città di un blocco lungo m
    void Individuo :: mutazione_inversione(){
    while(true){
        int max_m = _cities - 2;  // perché la prima città è fissa
        if (max_m < 2) return; // troppe poche città per fare l'inversione
        int m = static_cast<int>(_rnd->Rannyu(1, max_m)); //lunghezza blocco
        int start = static_cast<int>(_rnd->Rannyu(1, _cities - m)); //posizione di partenza
        if(m + start >= _cities) continue; //se la posizione del blocco e la sua lunghezza superano il numero delle città ripeti
        mat blocco = _individuo.rows(start, start + m - 1);
        blocco = flipud(blocco); // Inverte le righe
        _individuo.rows(start, start + m - 1) = blocco;//costruisco il nuovo individuo mutato

        break;
        }
        this->check();//controllo sulla prima città
    }

    //restituisce le riighe da a a b di individuo
mat Individuo::rows(int a, int b) {
    if (a > _cities - 1 || b > _cities - 1 || a < 0 || b < 0) {
        cout << "Le righe non esistono." << endl;
        return mat(); // ritorna matrice vuota
    }

    if (a > b) {
        cout << "Intervallo di righe non valido: a > b." << endl;
        return mat(); // ritorna matrice vuota
    }
    return _individuo.rows(a,b);

    }

    //restituisce la riga a di individuo
mat Individuo :: row(int a){
    if(a > _cities - 1 || a < 0){
        cout << "La riga non esiste." << endl;
        return mat(); // ritorna matrice vuota
    }

    return _individuo.row(a);

    }

    //restituisce il numero di città
int Individuo :: get_city_number(int i){
    return static_cast<int>(_individuo(i,0));

    }
//controllo che la città nella prima posizione sia quella etichettata con 1
void Individuo :: check(){
    if(_individuo(0,0)!=1.0)cout << "errore alle condizioni al contorno!" << endl;

    }
//restituisce l'individuo
mat Individuo :: get_individuo(){
    return _individuo;

}
//inizializza l'oggetto Random a cui il puntatore _rnd punta
void Individuo :: set_random(Random& rnd){ 
    _rnd = &rnd; 
}

