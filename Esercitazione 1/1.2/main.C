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

using namespace std;

//Prototipi:

double Err ( vector<double> & );				//Funzione di errore= dev stadard della media 
double Media ( vector<double> & );
double chi2 ( vector<int> &, int, int );
double chi1 ( double, double, double );


//Main 

int main () {

	Random R;															//Inizializzazione generatore casuale																		
	int seed[3];
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
				R.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else 
		cerr << "PROBLEM: Unable to open seed.in" <<endl;

	int M= 10000;
	vector <int> N = {1,2,10,100};
	
	ofstream out;												//ofstream per esportazione txt dei dati uniformi
	out.open("sim_uniforme_out.txt");
	ofstream out1;												//ofstream per esportazione txt dei dati esponenziali
	out1.open("sim_exp_out.txt");
	ofstream out2;												//ofstream per esportazione txt dei dati Lorentziani
	out2.open("sim_Lorentz_out.txt");
	
	
	for ( int j=0; j<M; j++) {
		for ( auto it:N ) {
			double a= 0.;
			double b= 0.;
			double c= 0.;
			for ( int i=0; i<it; i++ ) {			
				a+= R.uniforme();
				b+= R.exp( 1. );
				c+= R.Lorentz ( 1., 0. );
			}
			out <<fixed <<setprecision(6) << a/it << "\t";
			out1 <<fixed <<setprecision(6) << b/it << "\t";
			out2 <<fixed <<setprecision(6) << c/it << "\t";
 		}
 		out <<endl;
 		out1 <<endl;
 		out2 <<endl;
 	}	

	out.close();
	out1.close();
	out2.close();

	R.SaveSeed();
   
	return 0;
}


//Funzioni:

double Err ( vector<double> &A ) {
	
	if ( A.size()==0 || A.size()==1 )						//Controlli su valori problematici per operazioni di divisione [A.size() restituisce un intero non negativo]
		return 0;
	
	vector <double> A2;									//Vettore delle somme quadrate (cfr definizione di varianza)
	for ( auto it:A )
    	A2.push_back( pow(it,2));
	
	return sqrt( ( Media(A2) - pow(Media(A),2) ) / (A.size()-1) );						//Deviazione standard della media 
}

double Media ( vector <double> &X ) {

	return accumulate (X.begin(), X.end(), 0.)/X.size();

}

double chi2 ( vector<int> &X, int n, int M ) {

	if ( n==0 || M==0 )									//Controllo, anche se in verità non lo uso veramente perchè dovrei abortire il programma, invece do comunque un ritorno che può comuque avere un senso numerico ma non statistico
		return 0;
	
	double a= 0.;	
	for ( int i=0; i<M; i++ )
		a+= chi1( X[i], n, M );
		
	return a;

}

double chi1 ( double x, double n, double M ) {						//Arrivano int, ma li prendo double così non ho problemi con le divisioni

	return pow(x-n/M,2)/(n/M);

}


