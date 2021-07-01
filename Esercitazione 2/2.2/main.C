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
double Segno ( Random * );
void Passo ( vector<int> &, Random *, double );
vector<double> PassoContinuo ( vector<double> &, Random *, double );
vector<double> SphCar ( vector<double> & );
double Distanza1D ( double, double );
template <typename T> double Distanza ( vector<T>, vector<T> );


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

	vector<int> M;								//passi
	for ( int i=1; i<=100; i++)
		M.push_back (i);
	int N= 100;						//100 raggruppamenti su cui calcoliamo errore
	int L = 100;					//Dimensione dei raggruppamenti	
	double a= 1;							//passo reticolo 3D
	
	ofstream out;
	nomefile= "sim_discreto_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		for ( auto it:M ) {
			vector<double> Y;								//vector dove alloco la radice della distanza media di S ripetizioni dall'origine per ogni blocco di L elementi
			for ( int k=0; k<L; k++ ) {
				vector<int> X= {0,0,0};						//partenza della simulazione
				for ( int i=0; i<it; i++ )
					Passo ( X, R, a );
				Y.push_back ( pow(Distanza({0,0,0},X),2) );
			}
			out <<fixed <<setprecision(6) << Media(Y) << "\t";
		}
		out <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();
	nomefile= "sim_continuo_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		for ( auto it:M ) {
			vector<double> Y;								//vector dove alloco la radice della distanza media di S ripetizioni dall'origine per ogni blocco di L elementi
			for ( int k=0; k<L; k++ ) {
				vector<double> X= {0.,0.,0.};						//partenza della simulazione
				for ( int i=0; i<it; i++ )
					X= PassoContinuo ( X, R, a );
				Y.push_back ( pow(Distanza({0.,0.,0.},X),2) );
			}
			out <<fixed <<setprecision(6) << Media(Y) << "\t";
		}
		out <<endl;
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

double Segno ( Random* R ) {
	
	double x= R->uniforme(-1.,1.);
	
	if ( x==0 )
		return 1.;									//o si scarta x oppure arbitrariamente si sceglie +1 o -1, scelgo +1 così che sull'intervallo [-1,1) ci sia esattamente 50% di segni posotivi e negativi
	
	return x/fabs(x);								//sarebbe intero +1 o -1, ma lasciamo double prima che per qualche valore ci siano problemi 

}

void Passo ( vector<int> &A, Random *R, double a ) {

	int i= R->uniforme(0.,3.);
	A[i]+= a*Segno(R);
	
}

vector<double> PassoContinuo ( vector<double> &A, Random *R, double a ) {

		vector<double> B= { a,R->uniforme(0.,M_PI),R->uniforme(0.,2*M_PI) };
		vector<double> C= SphCar (B);
	
		return {A+C};
	}

vector<double> SphCar ( vector<double> &X ) {

	double x= X[0]*sin(X[1])*cos(X[2]);
	double y= X[0]*sin(X[1])*sin(X[2]);
	double z= X[0]*cos(X[1]);

	return {x,y,z};

}

double Distanza1D ( double x, double y ) {

	return fabs( x-y );
}

template <typename T> double Distanza ( vector<T> A, vector<T> B ) {

	if ( A.size() != B.size() )						//Controlli: non posso calcolare distanza tra due punti di due spazi dimensionali differenti
		return 0;
	
	double sum2= 0.;
	for ( unsigned int i=0; i<A.size(); i++ )
		sum2+= pow( Distanza1D ( A[i], B[i] ), 2);

	return sqrt( sum2 );
}


