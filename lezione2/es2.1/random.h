#ifndef __Random__
#define __Random__

// Questa classe contiene funzioni per generare numeri casuali usando l'algoritmo RANNYU
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Costruttore di default
  Random();
  // Distruttore
  ~Random();
  // Metodo per impostare il seme del generatore di numeri casuali
  void SetRandom(int * , int, int);
  // Metodo per salvare il seme su un file
  void SaveSeed();
  // Metodo per generare un numero casuale nell'intervallo [0,1)
  double Rannyu(void);
  // Metodo per generare un numero casuale nell'intervallo [min,max)
  double Rannyu(double min, double max);
  // Metodo per generare un numero casuale con distribuzione gaussiana
  double Gauss(double media, double sigma);
  // Metodo per generare un numero casuale con una distribuzione esponenziale
  double exp(double lambda);
  // Metodo per generare un numero casuale con una distribuzione di Cauchy-Lorentz
  double Cauchy(double mu, double gamma);
  //Metodo per generare un numero casuale secondo la distribuzione 2(x-1)
  double funzionesimil(); 
};

#endif // __Random__