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
double Prezzo ( double, double, double, double, double, Random * );
double Call ( double, double, double, double );
double Put ( double, double, double, double );


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
	int L = 100000;					//Dimensione dei raggruppamenti	
	
	double S0= 100.;
	double T= 1.;
	double K= 100.;
	double r= 0.1;
	double sigma= 0.25;


	/* Salto unico */

	ofstream out;
	nomefile= "sim_unico_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> C;								//vector dove alloco il guadagno delle Call del blocco da mediare
		vector<double> P; 								//vector dove alloco il guadagno delle Put del blocco da mediare
		for ( int k=0; k<L; k++ ) {
			double S= Prezzo (S0,0.,T,r,sigma,R);
			C.push_back ( Call(K,S,r,T) );
			S= Prezzo (S0,0.,T,r,sigma,R);						//Ridefiniamo perevitare che correlazioni anche se sono sati da analizzare separatamente
			P.push_back ( Put(K,S,r,T) );
		}
		out <<fixed <<setprecision(6) << Media(C) << "\t" << Media(P) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();


	/* Salti discreti */

	vector<double> Ti;								//passaggi discreti di tempo
	for ( int i=0; i<=100; i++)
		Ti.push_back (i/100.);
	
	nomefile= "sim_salti_out.txt";
	out.open( nomefile );

	for ( int j=0; j<N; j++ ) {
		vector<double> C;								//vector dove alloco il guadagno delle Call del blocco da mediare
		vector<double> P; 	
			for ( int k=0; k<L; k++ ) {
				double SC= S0;
				double SP= S0;
				for ( unsigned int i=1; i<Ti.size(); i++) {
					SC= Prezzo (SC,Ti[i-1],Ti[i],r,sigma,R);
					SP= Prezzo (SP,Ti[i-1],Ti[i],r,sigma,R);
				}
				C.push_back ( Call(K,SC,r,T) );
				P.push_back ( Put(K,SP,r,T) );	
			}		
		out <<fixed <<setprecision(6) << Media(C) << "\t" << Media(P) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
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

double Prezzo ( double S, double ti, double tf, double mean, double sigma, Random *R ) {

	double Z= R->Gauss(0.,1.); 
	
	return   S*exp( (mean-0.5*pow(sigma,2))*(tf-ti)+sigma*Z*sqrt(tf-ti) ); 

} 

double Call ( double K, double S, double mean, double T ) {

	return exp(-mean*T)*max(0., S-K);

}

double Put ( double K, double S, double mean, double T ) {

	return exp(-mean*T)*max(0., K-S);

}


