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
double Transizione ( double, double, Random * );
double phitrial ( double, double, double );
double Accettazione ( double, double, double, double );
double Hamiltoniana ( double, double, double );
double Potenziale ( double );

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
	int L = pow(10,6)/N;					//Dimensione dei raggruppamenti	t.c. si hanno 10^6 esecuzioni

	double mu= 0.82;
	double sigma= 0.64;
	double l= 2.75; 							//regola del 50% su segamnto di lato 2*l
	
	double x= 0.;													// punto di partenza della particella libera
	
	ofstream out;
	nomefile= "out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i valori del gruppo 
		int acc= 0; 
		for ( int k=0; k<L; k++ ) {
			double y= Transizione(x,l,R);
			if ( R->uniforme(0.,1.) <= Accettazione (x,y,mu,sigma) ) {
				x=y;
				acc++;
			}
			S.push_back ( Hamiltoniana(x,mu,sigma) );	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		cout << nomefile << ": " << j+1 << "%" 
			 << "con accettazione del Metropolis= " <<double(acc)/double(L) 
			 <<endl;
	}
	cout <<endl;
	
	out.close();
				
	R->SaveSeed();
   
	return 0;
}


//Funzioni:

double Media ( vector <double> &A ) {

	return accumulate ( A.begin(), A.end(), 0. )/A.size();

}

double Transizione ( double a, double l, Random *R ) {

	return a + R->uniforme( (-1.)*l,l );

}

double phitrial ( double x, double mu, double sigma ) {

	return exp( (-1.)*pow(x-mu,2)/(2*pow(sigma,2)) ) + exp( (-1.)*pow(x+mu,2)/(2*pow(sigma,2)) );

}


double Accettazione ( double xold, double xnew, double mu, double sigma ) {
	
	return min( 1., pow(phitrial(xnew,mu,sigma),2)/pow(phitrial(xold,mu,sigma),2) );

}

double Hamiltoniana ( double x, double mu, double sigma ) {

	double phi= phitrial(x,mu,sigma);
	double Kin= 1./(2*pow(sigma,2)) * ( phi - pow((x-mu)/sigma,2)*exp( (-1.)*pow(x-mu,2)/(2*pow(sigma,2)) ) - pow((x+mu)/sigma,2)*exp( (-1.)*pow(x+mu,2)/(2*pow(sigma,2)) ) );
	return ( Kin + Potenziale(x)*phi )/phi;

}

double Potenziale ( double x ) {											//V(x)= x^4 - 5/2 x^2

	 return pow(x,4) - 5./2.*pow(x,2);

}




