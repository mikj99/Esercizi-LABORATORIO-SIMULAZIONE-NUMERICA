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

	int N= 100;														//100 raggruppamenti
	int L= 1000;														//Dimensione dei raggruppamenti	
	
	vector<double> X;											//Variabili necessarie allo svolgimento
	vector<double> X2;
	vector<double> err_prog;
	vector<double> sigma2;
	
	ofstream out;												//ofstream per esportazione txt
	out.open("sim_out.txt");
	
	for ( int j=0; j<N; j++ ) {									//Generazione 100 raggruppamenti come media di L numeri casuali
		double sum=0.;
		double sum2=0.;
		for ( int i=0; i<L; i++ ){			
			double a= R.Rannyu();
			sum+= a;										//Generazione numeri casuali			
 			sum2+=pow(a-0.5,2);
 		}
		X.push_back( sum/L );												//Vettore contenete le misure
		X2.push_back( sum2/L );
		err_prog.push_back( Err(X) );															//Vettore degli errori progessivi sulle misure contenute nel vettore X
		sigma2.push_back ( Err(X2) );
		out <<fixed <<setprecision(6) << Media(X) << "\t" << err_prog.back() << "\t" << Media(X2) << "\t" << sigma2.back() <<endl;										//txt in output
	}
	
	out.close();
	out.open("chi2.txt");
	
	int M= 100;
	int n= 10000;
	
	for ( int j=0; j<N; j++ ) {
		vector<int> span (M);
		for ( int i=0; i<n; i++) 
			span[ int( R.Rannyu()*M ) ]++;
		out <<fixed <<setprecision(2) << chi2(span, n, M) <<endl;
	}

	out.close();

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


