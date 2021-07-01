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

	
	int N= 5*pow(10,4);
	double mu, sigma, E=1000.;
	double l= 2.75;					//usiamo un l che che permetta una accettazione del Metropolis del 50% per il minimo trovato a mano dell'esercizio 8.1, 
							//ci muoveremo in un suo intorno, pertanto la accettazione dovrebbe rimanere intono al 50%  	
	
	for ( int i=0; i<10000; i++ ) {
		double x= 0.;
		double mu1= R->uniforme(.77,.87);
		double sigma1= R->uniforme(.57,.67);	
		vector <double> A;						//vector dove alloco le stime di energia per ogni sim
		for ( int k=0; k<N; k++ ) {
			double y= Transizione(x,l,R);
			if ( R->uniforme(0.,1.) <= Accettazione (x,y,mu1,sigma1) )
				x=y;
			A.push_back ( Hamiltoniana(x,mu1,sigma1) );	
		}		
		if ( Media(A) < E ) {
			mu= mu1;
			sigma= sigma1;
			E= Media(A);
			cout << "Migliori parametri trovati ad ora:" <<endl;
			cout <<fixed <<setprecision(4) << "- mu    = " <<mu <<endl;
			cout <<fixed <<setprecision(4) << "- sigma = " <<sigma <<endl;
			cout <<fixed <<setprecision(6) << "- E0    = " <<E <<endl;
			cout << "-----------------------------------------" <<endl <<endl;
		}
	}
	
	//simulazione con i migliori parametri:

	ofstream out;
	nomefile= "parametri.txt";
	out.open( nomefile );
	cout << "Stampo nel file " <<nomefile << "i migliori parametri ottenuti secondo l'ordine mu|sigma|E0." <<endl <<endl;
	out <<fixed <<setprecision(6) <<mu << "\t" <<sigma << "\t" <<E <<endl;
	out.close(); 
	

	N= 100;						//100 raggruppamenti su cui calcoliamo errore
	int L = pow(10,6)/N;					//Dimensione dei raggruppamenti	t.c. si hanno 10^6 esecuzioni
	
	double x= 0.;													// punto di partenza della particella libera
	
	nomefile= "out.txt";
	out.open( nomefile );
	ofstream out1;
	string nomefile1= "funzioned'onda.txt";
	out1.open( nomefile1 );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i valori del gruppo 
		for ( int k=0; k<L; k++ ) {
			double y= Transizione(x,l,R);
			if ( R->uniforme(0.,1.) <= Accettazione (x,y,mu,sigma) )
				x=y;
			S.push_back ( Hamiltoniana(x,mu,sigma) );
			out1 <<fixed <<setprecision(6) <<x <<endl;	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
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




