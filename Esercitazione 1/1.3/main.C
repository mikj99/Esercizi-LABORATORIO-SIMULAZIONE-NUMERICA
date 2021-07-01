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
double Distanza1D ( double, double );
double Distanza ( vector<double> &, vector<double> & );
bool Intersezione ( double, vector<double> &, vector<double> & );
vector<double> PuntoMedio ( vector<double> &, vector<double> & );
vector<double> Estremo2 ( vector<double> &, double, Random * );
double seno ( double, Random * );

//Main 

int main () {

	Random* R= new Random ();			//Inizializzazione generatore casuale																		
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
				R->SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else 
		cerr << "PROBLEM: Unable to open seed.in" <<endl;

	double d= 5.;				   //aste a distanza 2 lingo y
	double rho= 2.5;				//lunghezza barrette lanciate
	double min= 10.;					//limiti dello spazio --> (x,y) = ( [10,30], [10,30] ) mi sposto dall'origine sennò int(-0.5)=0 e non funziona Intersezione ( double, vector<double>, vector<double> )
	double max= 30.;
	double M=pow(10,3);					//Lanci per stima singolo pi
	
	
/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	
	/* METODO PUNTO MEDIO (hit rate ~ 100.4% su 10^6 ripetizioni ==> leggera sovrastima che tende adaccentuarsi nel numero di tentativi)
		int n_hit= 0;
		for ( int i=0; i<M; i++ ) {
			vector<double> A;
			vector<double> B;
			for ( int k=0; k<2; k++ ) {				//R2, spazio 2 dimensionale, genero Xa, Xb, Ya, Yb, così che X e T dello stesso punto non siano menneno consecutivi
				A.push_back( R->uniforme(min, max) );
				B.push_back( R->uniforme(min, max) );	
			}
			vector<double> M= PuntoMedio (A,B);
			vector<double> C;
			vector<double> D;
			for ( unsigned int i=0; i<A.size(); i++ ) {
				C.push_back( M[i]- (B[i]-A[i])*rho/(2*Distanza(A,B)) );
				D.push_back( M[i]+ (B[i]-A[i])*rho/(2*Distanza(A,B)) );
			}
				
			if ( Intersezione(d,C,D) )
				n_hit++;
		}	
		cout << "pi= " <<(2*rho*M)/(n_hit*d) << "\tn_hit= " <<n_hit << "\th_hit_teo= " << (2*rho*M)/(M_PI*d) <<endl;
	*/
	
	/* METODO RETTA SEMPLICE (hit rate ~ 93.3% su 10^6 ripetizioni ==> sottostima che diminuisce nel numero di tentativi)
		int n_hit= 0;
		for ( int i=0; i<M; i++ ) {		
			vector<double> A;
			vector<double> B;
			for ( int k=0; k<2; k++ ) {				//R2, spazio 2 dimensionale, genero Xa, Xb, Ya, Yb, così che X e T dello stesso punto non siano menneno consecutivi
				A.push_back( R->uniforme(min, max) );
				B.push_back( R->uniforme(min, max) );	
			}
			vector<double> C;
			for ( unsigned int i=0; i<A.size(); i++ )
				C.push_back( A[i]+ (B[i]-A[i])*rho/(Distanza(A,B)) );
				
			if ( Intersezione(d,A,C) )
					n_hit++;
		}
		cout << "pi= " <<(2*rho*M)/(n_hit*d) << "\tn_hit= " <<n_hit << "\th_hit_teo= " << (2*rho*M)/(M_PI*d) <<endl;
	*/
	
	/* METODO PROIEZIONE COSENO (hit rate ~ 99.88-100.06% su 10^6 ripetizioni ==> stima oscillante) 
		int n_hit= 0;
		for ( int i=0; i<M; i++ ) {
			vector<double> A;
			for ( int k=0; k<2; k++ )				//R2, spazio 2 dimensionale
				A.push_back( R->uniforme( min, max ) );	
			vector<double> B= Estremo2 ( A,rho,R );
				
			if ( Intersezione( d,A,B ) )
				n_hit++;
		}	
		cout <<fixed <<setprecision(6) << "pi= " <<(2*rho*M)/(n_hit*d) << "\tn_hit= " <<n_hit << "\th_hit_teo= " << (2*rho*M)/(M_PI*d) <<endl;
	*/
	
/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	
	
	int N= 100;						//100 raggruppamenti su cui calcoliamo errore
	vector <int> L = {1,10,100,1000};					//Dimensione dei raggruppamenti	
	
	ofstream out;												//ofstream per esportazione txt
	out.open("sim_out.txt");
	
	for ( int j=0; j<N; j++ ) { 
		for ( auto it:L ) {
			vector<double> Y;
			for ( int k=0; k<it; k++ ) {
				int n_hit= 0;
				for ( int i=0; i<M; i++ ) {
					vector<double> A;
					for ( int k=0; k<2; k++ )						//R2, spazio 2 dimensionale
						A.push_back( R->uniforme( min, max ) );	
					
					vector<double> B= Estremo2 ( A,rho,R );
					if ( Intersezione( d,A,B ) )
						n_hit++;
				}
				Y.push_back( (2*rho*M)/(n_hit*d) );
			}
			out <<fixed <<setprecision(6) << Media(Y) << "\t";				//output ci sono solo i valori medi dei gruppi, nel notebook si farà l'andamento delle somme e la varianza
		}
		out <<endl;
		cout << j+1 << "%" <<endl;
	}
	
	out.close();
	
	R->SaveSeed();
   
	return 0;
}


//Funzioni:

double Media ( vector <double> &X ) {

	return accumulate ( X.begin(), X.end(), 0. )/X.size();

}

double Distanza1D ( double x, double y ) {

	return fabs( x-y );
}

double Distanza ( vector<double> &A, vector<double> &B ) {

	if ( A.size() != B.size() )						//Controlli: non posso calcolare distanza tra due punti di due spazi dimensionali differenti
		return 0;
	
	double sum2= 0.;
	for ( unsigned int i=0; i<A.size(); i++ )
		sum2+= pow( Distanza1D ( A[i], B[i] ), 2);

	return sqrt( sum2 );
}

bool Intersezione ( double d, vector<double> &A, vector<double> &B ) {					//Funziona solo per il nostro esperimento: A e B in R2, d lungo y ([1])

	if ( int(A[1]/d) == int(B[1]/d) )
		return false;
			
	return true;
	
}

vector<double> PuntoMedio ( vector<double> &A, vector<double> &B ) {

	if ( A.size() != B.size() )						//Controlli: non posso cercare punto medio di un segmento i cui gli estremi sono appartenti a spazi dimensionali diversi
		return {0.};
	
	vector<double> M;
	for ( unsigned int i=0; i<A.size(); i++ )
		M.push_back( A[i]+ (B[i]-A[i])/2 );

	return M;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

	/* Utilizzo theta  [0,2pi]

	vector<double> Estremo2 ( vector<double> &X, double L, Random *R ) {

		vector<double> Y;
		double theta= R->uniforme(0.,2*M_PI);
		Y.push_back( X[0] + L*cos(theta) );
		Y.push_back( X[1] + L*seno(cos(theta),R) );
 
		return Y;
	}
	*/

/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


vector<double> Estremo2 ( vector<double> &X, double L, Random *R ) {						//in R2, caso specifico del problema

	vector<double> Y;
	double cos= R->coseno();
	Y.push_back( X[0] + L*cos );
	Y.push_back( X[1] + L*seno(cos,R) );
 
	return Y;
}

double seno ( double cos, Random *R ) {
	
	double segno= R->uniforme(-1,1);
	return ( segno*sqrt(1-pow(cos,2))/fabs(segno) );

}



