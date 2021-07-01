#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "overloading.h"

using namespace std;

// Prototipi:

void LoadFile ( vector<vector<double>> &, string );
void Quadrato ( vector<vector<double>> &, int, double, Random * );
void Circonferenza ( vector<vector<double>> &, int, double, Random * );
void Move ( vector<int> &, Random * );
void Accetta ( vector<vector<double>>, vector<int> &, vector<int> &, Random *, double );
double FLoss ( vector<vector<double>>, vector<int> );
double distanza ( vector<double>, vector<double> );
void BestPath ( vector<vector<double>>, vector<int>, string );

// Main

int main () {

	Random* R= new Random ();								//set random number generator
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else {
		cerr << "PROBLEMA: Incapace di aprire Primes" << endl;
		exit(1);
	}
	Primes.close();

	ifstream input( "seed.out" );
	if ( input.is_open() ){
		while ( !input.eof() ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			R->SetRandom(seed,p1,p2);
			}
		input.close();
	} else {
		cerr << "PROBLEMA: Incapace di aprire seed.out; provo con seed.in" <<endl;
		input.close();
		input.open( "seed.in" );
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
		} else {
			cerr << "PROBLEMA: Incapace di aprire seed.in" <<endl;
			exit(1);
		}
	}

	
	cout <<endl << "Scegliere in quale spazio estarre le città da visitare" <<endl;
	cout << "[quadrato/circonferenza/file esterno]" <<endl;
	string spazio;
	cin >>spazio;
	if ( spazio!="quadrato" && spazio!="circonferenza" && spazio!="file esterno" ) {
		cout << "Conformazione dello spazio non conforme, rieseguire il programma e seguire attentamente le istruzioni" <<endl;
		exit(1);
	}
	
	int n;
	if ( spazio=="quadrato" || spazio=="circonferenza" ) {
		cout <<endl << "Scegliere il numero (intero) delle città da visitare nello spazio scelto, essendo questo non proveniente da un file esterno" <<endl;
		cin >>n;
		cout <<endl;
	} 

	string nome;
	if ( spazio=="file esterno" ) {
		cout <<endl << "Immetere il nome del file esterno (solo 2-dimensionale) da cui prendere la configurazione delle città da visitare." <<endl;
		cout << "(Verranno utilizzate tutte le città riporate, qualora il file non sia 2-dimensionale il programma produrrà risulatati non corretti e controllati)" <<endl; 
		cin >>nome;
		cout <<endl;
	}
	
	double beta= 1.;													//beta= 1/T
	cout << "Scegliere per quante generazioni (numero intero) eseguire l'evoluzione (almeno 200)" <<endl;				//cambio T ad ogni generazione
	int generazioni;
	cin >>generazioni;
	cout <<endl;
	
	vector<vector<double>> X;						//vettore delle posizioni delle città
	if ( spazio=="quadrato" )
		Quadrato( X, n, 10., R );
	if ( spazio=="circonferenza" )
		Circonferenza( X, n, 1., R );
	if ( spazio=="file esterno" )
		LoadFile( X, nome );

	vector<int> Cromo;
	for ( unsigned int i=0; i<X.size(); i++ )
		Cromo.push_back(i);
	vector<int> Cromo1= Cromo;
	int pop= 200000;
	
	ofstream fl;
	fl.open( "FLoss_" + spazio + ".txt" );
	cout << "Stampo nel file FLoss_" <<spazio << ".txt alcuni dati relativi alla funzione obiettivo; in paticolare:\ngenerazione | temperatura | beta | funzione obiettivo "<<endl <<endl;
	
	cout << "Eseguo per " <<generazioni << " generazioni:" <<endl;
	for ( int i=1; i<=generazioni; i++ ) {
		
		for ( int j=0; j<pop; j++ ) {
		
			Move ( Cromo1,R );
			Accetta ( X,Cromo,Cromo1,R,beta ); 
		
		}
		
		fl <<i << "\t" << 1./beta << "\t" <<beta << "\t" <<FLoss(X,Cromo) <<endl;
		
		if ( i%int(double(generazioni)/100.) == 0 )
			cout <<double(i)/double(generazioni)*100 << "%\t" <<flush;
		if ( i%int(double(generazioni)/10.) == 0 )
			cout <<endl;

		beta*= 1.03;
		pop= 200000./beta + 0.5;
	
	}
	cout<<endl;
	
	fl.close();
	BestPath( X,Cromo, "bestpath_" + spazio + ".txt" );
	
	R->SaveSeed();
	   
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// FUNZIONI:

void LoadFile ( vector<vector<double>> &X, string name ) { 
	
	ifstream Input( name );
	if ( Input.is_open() ) {
		while ( !Input.eof() ){
			double x,y;
			Input >> x >> y;
			X.push_back( {x,y} );
		}
	} else {
		cerr << "PROBLEMA: impossibile caricare la configurazione iniale dal file fornito " <<name << endl;
		exit(1);
	}
	Input.close();
	
}

void Quadrato ( vector<vector<double>> &X, int n, double l, Random *R ) {
	
	for ( int i=0; i<n; i++) {
		double x= R->uniforme(0.,l);
		double y= R->uniforme(0.,l);
		X.push_back( {x,y} );
	}

}

void Circonferenza ( vector<vector<double>> &X, int n, double r, Random *R ) {
	
	for ( int i=0; i<n; i++) {
		double theta= R->uniforme(0.,2.*M_PI);
		X.push_back( {r*cos(theta),r*sin(theta)} );
	}

}

void Move ( vector<int> &A, Random *R ) {				//permutazione semplice
		
	int a= R->uniforme(0.,4.);
	
	if ( a == 0 ) {
		int i= R->uniforme(1.,A.size()), j= R->uniforme(1.,A.size());	
		swap( A[i], A[j] );
	}

	if ( a == 1 ) {
		vector<int> Y;
		for ( int i=0; i<3; i++ )
			Y.push_back ( R->uniforme(1.,A.size()) );
		sort( Y.begin(), Y.end() );
		
		rotate( A.begin()+Y[0], A.begin()+Y[1], A.begin()+Y[2] );
	}

	if ( a == 2 ) {
		int i= R->uniforme(1.,A.size()), j= R->uniforme(1.,A.size());
		if ( i>j )
			swap( i,j );
		int m= R->uniforme(1.,A.size()-j);
		for ( int k=0; k<m; k++ )
			swap( A[i+k],A[j+k] );
	}
	
	if ( a == 3 ) {
		int i= R->uniforme(1.,A.size()), j= R->uniforme(1.,A.size());	
		if ( i>j )
			swap( i,j );
		while ( i<j ) {
			swap( A[i],A[j] );
			i++;
			j--;
		}
	}
	
}

void Accetta ( vector<vector<double>> X, vector<int> &A, vector<int> &B, Random *R, double beta ) {

	double delta= FLoss(X,B)-FLoss(X,A); 
	
	if ( delta < 0 )
		A= B;
	else
		if ( R->uniforme(0.,1.) < exp( (-1.)*beta*delta ) )
			A= B;
		else 
			B= A;

}

double FLoss ( vector<vector<double>> X, vector<int> Y ) {

	double d= 0.;
	for ( unsigned int i=0; i<Y.size()-1; i++ )
		d+= distanza( X[Y[i]], X[Y[i+1]] );
	d+= distanza ( X[Y[0]], X[Y[Y.size()-1]] );

	return d;
	
}

double distanza ( vector<double> B, vector<double> C ) {

	//Controlli
	if ( B.size() != C.size() ){
		cout << "Richiesta una distanza tra oggeti in spazi dimensionalmente diversi!\nErrore!" <<endl;
		exit(1);
	}
	if ( B.size() == 0 || C.size() == 0 ){
		cout << "Richiesta una distanza tra oggeti in spazi dimensionalmente nulli!\nErrore!" <<endl;
		exit(1);
	}

	double s2= 0.;
	for ( unsigned int i=0; i<B.size(); i++ )
		s2+= pow(B[i]-C[i],2);
	
	return sqrt( s2 );

}

void BestPath ( vector<vector<double>> X, vector<int> Y, string name ) {

	cout << "Stampo sul file " <<name << " le posizioni degli oggetti allocati nell'ultimo cromosoma ottenuto" <<endl;
	
	ofstream out;
	out.open( name );
	if ( out.is_open() ) {
		for ( auto ity : Y ) {
			out << ity << "\t";
			for ( auto itx : X[ity] )
				out <<fixed <<setprecision(6) <<itx << "\t";
			out <<endl;
		}
	} 
	else
		cerr << "PROBLEMA: Impossibile aprire posizioni.txt" << endl;
	out.close();
	
}


