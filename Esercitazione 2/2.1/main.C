#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include "random.h"
#include "overloading.h"

using namespace std;

//Prototipi:

double Media ( vector<double> & );
double AccRej ( double, double, double, Random * );



//Main 

int main () {

	Random* R= new Random ();										//Inizializzazione generatore casuale																		
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	string nomefile= "seed.out";
	ifstream input( nomefile );
	if ( input.is_open() ){
		while ( !input.eof() ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			R->SetRandom(seed,p1,p2);
			}
		input.close();
	} else {
		cerr << "PROBLEM: Unable to open seed.out; i'm trying seed.in" <<endl;
		input.close();
		nomefile= "seed.in";
		input.open( nomefile );
		string property;
		if ( input.is_open() ){
			while ( !input.eof() ){
				input >> property;
				if( property == "RANDOMSEED" ){
					input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					R->SetRandom(seed,p1,p2);
				}
			}
			input.close();
		} else 
		cerr << "PROBLEM: Unable to open seed.in" <<endl;
	}

	int N= 100;						//100 raggruppamenti su cui calcoliamo errore
	vector <int> L = {1,10,100,1000};					//Dimensione dei raggruppamenti	
	
	ofstream out;												//ofstream per esportazione txt
	out.open("sim_uniforme_out.txt");
	
	for ( int j=0; j<N; j++ ) { 						//p(x) uniforme [0,1)
		for ( auto it:L ) {
			vector<double> Y;
			for ( int k=0; k<it; k++ )
				Y.push_back( cos(M_PI*R->uniforme(0.,1.)/2) );
			out <<fixed <<setprecision(6) << M_PI*Media(Y)/2 << "\t";
		}
		out <<endl;
		cout << "sim_uniforme_out.txt  " << j+1 << "%" <<endl;
	}
	cout <<endl;
	out.close();
	
	out.open("sim_p(x)_out.txt");
	
	for ( int j=0; j<N; j++ ) { 								//p(x)= 1-x^2  g(x)= 2cos(x)/(3*(1-x^2)) con pmax=1
		for ( auto it:L ) {
			vector<double> Y;
			for ( int k=0; k<it; k++ ) {
				double x= AccRej(0.,1.,1.,R);
				Y.push_back( 2.*cos(M_PI*x/2.)/(3.*(1.-pow(x,2.))) );		
				}
			out <<fixed <<setprecision(6) << M_PI*Media(Y)/2 << "\t";				//output ci sono solo i valori medi dei gruppi, nel notebook si farà l'andamento delle somme e la varianza
		}
		out <<endl;
		cout << "sim_p(x)_out.txt  " << j+1 << "%" <<endl;
	}
	cout <<endl;
	out.close();
	
	R->SaveSeed();
   
	return 0;
}


//Funzioni:

double Media ( vector <double> &X ) {

	return accumulate ( X.begin(), X.end(), 0. )/X.size();

}

double AccRej ( double min, double max, double pmax, Random *R ) {
	  
	  double x;
	  double y;
	  
	  do {
	  x= R->uniforme(min,max);
	  y= R->uniforme(0.,pmax);
	  } while ( y> (1-pow(x,2)) );
	  
	  return x;

}


