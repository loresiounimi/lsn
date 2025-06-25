#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
  string isto[3]= {"standardisto.txt", "expisto.txt", "cauchyisto.txt"};
  ofstream istowrite[3];
  for(int i=0;i<3;i++){
    istowrite[i].open(isto[i]);
  if(!istowrite[i]){
    cout << "Error: Unable to open the file" << endl;
    return 1;
    } 
    } 
  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
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
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  int nmax= 10000;
  double dep=0.0;
  vector <int> s{1,2,10,100};
  vector <vector<vector<double>>> v(3);
  for(int y=0; y<3; y++){
    v[y].resize(s.size());
  for(int j=0; j<(int)s.size(); j++){
  for(int i=0; i<nmax; i++){
    for(int h=0; h<s[j]; h++){
      if(y==0){
        dep+=rnd.Rannyu();
        }else if(y==1){
        dep+=rnd.exp(1.0);
        }else if(y==2){
        dep+=rnd.Cauchy(0.0,1.0);
        }else{
        break;
        }
      }
    v[y][j].push_back(dep/s[j]);
    dep=0.0;
    }
    }
  for(int i=0; i<nmax; i++){
    istowrite[y] << v[y][0][i] << " " << v[y][1][i] << " " << v[y][2][i] << " " << v[y][3][i] << endl;
    }
    //v[y].clear();
    }
  return 0;
}