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
vector<double> SphCar ( vector<double> );
double raggio ( vector<double> );
vector<double> Transizione ( vector<double>, double, Random * );
vector<double> TransizioneGauss ( vector<double>, double, Random * );
double Accettazionephi100 ( vector<double>, vector<double> );
double Accettazionephi210 ( vector<double>, vector<double> );


//Variabili genarali
double a= 0.0529*pow(10,-9);								//a0

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

	// Stato fondamentale T(x|y) uniforme
	vector<double> X= {1.,0.,0.};													// vettore di partenza con raggio in unità di a0
	double step= 1.225;																//regola del 50% su un cubo di lato 2*step
	
	ofstream out;
	nomefile= "fondamentale_uniforme_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i raggi nel gruppo
		for ( int k=0; k<L; k++ ) {
			vector<double> Y= Transizione(X,step,R);
			if ( R->uniforme(0.,1.) <= Accettazionephi100 (Y,X) )
				X=Y;
			S.push_back ( raggio(X) );	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		if ( j==0 ) {
			ofstream out1 ( "covarianza_out.txt" );
			for ( auto it:S ) 
				out1 <<fixed <<setprecision(6) <<it <<endl;
			out1.close();
		}
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();
	
	// Primo stato eccitato  T(x|y) uniforme
	X= {1.,0.,0.};
	step= 2.98;
	
	nomefile= "1eccitato_uniforme_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i raggi nel gruppo
		for ( int k=0; k<L; k++ ) {
			vector<double> Y= Transizione(X,step,R);
			if ( R->uniforme(0.,1.) <= Accettazionephi210 (Y,X) )
				X=Y;
			S.push_back ( raggio(X) );	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();
	
	// Stato fondamentale T(x|y) gaussiano
	X= {1.,0.,0.};													// vettore di partenza con raggio in unità di a0
	double sigma= 0.76;																//regola del 50% su un cubo di lato 2*step
	
	nomefile= "fondamentale_gauss_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i raggi nel gruppo
		for ( int k=0; k<L; k++ ) {
			vector<double> Y= TransizioneGauss(X,sigma,R);
			if ( R->uniforme(0.,1.) <= Accettazionephi100 (Y,X) )
				X=Y;
			S.push_back ( raggio(X) );	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();
	
	// Primo stato eccitato  T(x|y) gaussiano
	X= {1.,0.,0.};
	sigma= 1.88;
	
	nomefile= "1eccitato_gauss_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i raggi nel gruppo
		for ( int k=0; k<L; k++ ) {
			vector<double> Y= TransizioneGauss(X,sigma,R);
			if ( R->uniforme(0.,1.) <= Accettazionephi210 (Y,X) )
				X=Y;
			S.push_back ( raggio(X) );	
		}		
		out <<fixed <<setprecision(6) << Media(S) <<endl;
		cout << nomefile << ": " << j+1 << "%" <<endl;
	}
	cout <<endl;
	
	out.close();
		
	// Fondamentale  T(x|y) gaussiano con partenza lontano dall'origine
	X= {10.,0.,0.};
	L= 3*pow(10,4)/N;
	sigma= .76;
	
	nomefile= "fondamentale_gauss_lontano_out.txt";
	out.open( nomefile );
	
	for ( int j=0; j<N; j++ ) {
		vector<double> S;														//vector dove alloco i raggi nel gruppo
		for ( int k=0; k<L; k++ ) {
			vector<double> Y= TransizioneGauss(X,sigma,R);
			if ( R->uniforme(0.,1.) <= Accettazionephi100 (Y,X) )
				X=Y;
			S.push_back ( raggio(X) );	
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

double Media ( vector <double> &X ) {

	return accumulate ( X.begin(), X.end(), 0. )/X.size();

}

vector<double> SphCar ( vector<double> X ) {

	double x= X[0]*sin(X[1])*cos(X[2]);
	double y= X[0]*sin(X[1])*sin(X[2]);
	double z= X[0]*cos(X[1]);

	return {x,y,z};

}

double raggio ( vector<double> A ) {
	double s= 0.;
	for ( unsigned int i=0; i<A.size(); i++)
		s+= pow(A[i],2);

	return sqrt( s );

}

vector<double> Transizione ( vector<double> A, double step, Random *R ) {

	vector<double> B= { R->uniforme(-step,step), R->uniforme(-step,step), R->uniforme(-step,step) };
	return A+B;

}

vector<double> TransizioneGauss ( vector<double> A, double sigma, Random *R ) {

	return { R->Gauss(A[0],sigma), R->Gauss(A[1],sigma), R->Gauss(A[2],sigma) };

}


double Accettazionephi100 ( vector<double> A, vector<double> B ) {
	
	return min( 1., exp( 2*(raggio(B)-raggio(A)) ) );														// raggio in unità di a0

}

double Accettazionephi210 ( vector<double> A, vector<double> B ) {

	return min( 1., pow(A[2]/B[2],2)*exp( raggio(B)-raggio(A) ) );

}


