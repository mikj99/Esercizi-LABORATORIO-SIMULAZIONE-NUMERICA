#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include "random.h"
#include "genetic.h"

using namespace std;

//Protected:

double Genetic :: Media ( vector <double> A ) {

	return accumulate ( A.begin(), A.end(), 0. )/A.size();

}

double Genetic :: distanza ( vector<double> B, vector<double> C ) {

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

double Genetic :: FLoss ( vector<int> Y ) {

	double d= 0.;
	for ( int i=0; i<ngeni-1; i++ )
		d+= distanza( X[Y[i]], X[Y[i+1]] );
	d+= distanza ( X[Y[0]], X[Y[ngeni-1]] );

	return d;
	
} 

bool Genetic :: minFLoss ( vector<int> B, vector<int> C ) {

	return ( FLoss(B) < FLoss(C) );
	
}

void Genetic :: Popolare () {

	for ( unsigned int i=0; i<X.size(); i++ )
		A.push_back(i);
	ngeni= A.size();
	npop= pow(A.size(),2);
		
	for ( int i=0; i<npop; i++ ) {
		vector<int> B= A, C{0};
		while ( B.size()>1 ) {
			int k= R->uniforme(1.,B.size());
			C.push_back( B[k] );
			B.erase( B.begin()+k );
		}	
		P.push_back ( C );
	}
	ngen= 0;
	Ordine ();
	
}

void Genetic :: Mutare1 () {				//permutazione semplice
	
	int muta= R->uniforme(0.,4.);
	//muta= 0 ---> non mutano ne F1 ne F2
	//muta= 1 ---> muta solo F1
	//muta= 2 ---> muta solo F2
	//muta= 3 ---> mutano F1 ed F2
	//F1 ed F2 mutano il 50% delle volte che viene chiamato Mutare1(), quindi il 10% delle volte che viene chiamato Mutare()
	
	if ( muta == 0 )
		return;
		
	if ( muta == 1 || muta == 3 ) {	
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);	
		swap( F1[i],F1[j] );
	}
	
	if ( muta == 2 || muta == 3 ) {
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);
		swap ( F2[i],F2[j] );
	}
		
}

void Genetic :: Mutare2 () {				//sposto m elementi a partire dalla posizione i di +n

	int muta= R->uniforme(0.,4.);
	//muta= 0 ---> non mutano ne F1 ne F2
	//muta= 1 ---> muta solo F1
	//muta= 2 ---> muta solo F2
	//muta= 3 ---> mutano F1 ed F2
	//F1 ed F2 mutano il 50% delle volte che viene chiamato Mutare2(), quindi il 10% delle volte che viene chiamato Mutare()
	
	if ( muta == 0 )
		return;
		
	if ( muta == 1 || muta == 3 ) {	
		vector<int> Y;
		for ( int i=0; i<3; i++ )
			Y.push_back ( R->uniforme(1.,ngeni) );
		sort( Y.begin(), Y.end() );
		
		rotate( F1.begin()+Y[0], F1.begin()+Y[1], F1.begin()+Y[2] );
	}
	
	if ( muta == 2 || muta == 3 ) {
		vector<int> Y;
		for ( int i=0; i<3; i++ )
			Y.push_back ( R->uniforme(1.,ngeni) );
		sort( Y.begin(), Y.end() );
		
		rotate( F2.begin()+Y[0], F2.begin()+Y[1], F2.begin()+Y[2] );
	}
		
}

void Genetic :: Mutare3 () {					//permutazione di blocco contiguo di m elementi

	int muta= R->uniforme(0.,4.);
	//muta= 0 ---> non mutano ne F1 ne F2
	//muta= 1 ---> muta solo F1
	//muta= 2 ---> muta solo F2
	//muta= 3 ---> mutano F1 ed F2
	//F1 ed F2 mutano il 50% delle volte che viene chiamato Mutare3(), quindi il 10% delle volte che viene chiamato Mutare()
	
	if ( muta == 0 )
		return;
		
	if ( muta == 1 || muta == 3 ) {	
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);
		if ( i>j )
			swap( i,j );
		int m= R->uniforme(1.,ngeni-j);
		for ( int k=0; k<m; k++ )
			swap( F1[i+k],F1[j+k] );
	}
	
	if ( muta == 2 || muta == 3 ) {
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);
		if ( i>j )
			swap( i,j );
		int m= R->uniforme(1.,ngeni-j);	
		for ( int k=0; k<m; k++ )
			swap( F2[i+k],F2[j+k] );
	}
	
}

void Genetic :: Mutare4 () {					//inversione dei geni in [i,j] con i<=j

	int muta= R->uniforme(0.,4.);
	//muta= 0 ---> non mutano ne F1 ne F2
	//muta= 1 ---> muta solo F1
	//muta= 2 ---> muta solo F2
	//muta= 3 ---> mutano F1 ed F2
	//F1 ed F2 mutano il 50% delle volte che viene chiamato Mutare4(), quindi il 10% delle volte che viene chiamato Mutare()
	
	if ( muta == 0 )
		return;
		
	if ( muta == 1 || muta == 3 ) {	
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);	
		if ( i>j )
			swap( i,j );
		while ( i<j ) {
			swap( F1[i],F1[j] );
			i++;
			j--;
		}
	}
	
	if ( muta == 2 || muta == 3 ) {
		int i= R->uniforme(1.,ngeni), j= R->uniforme(1.,ngeni);	
		if ( i>j )
			swap( i,j );
		while ( i<j ) {
			swap( F2[i],F2[j] );
			i++;
			j--;
		}
	}
	
}

vector<int> Genetic :: Soluzione () {

	Ordine();
	return P[0];

}

/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//Public:

Genetic :: Genetic(){
	
	//Set generatore di numeri casuali
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

}

Genetic :: ~Genetic() {

	R->SaveSeed();

}

void Genetic :: SetGenetic ( string name ) { 

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
	
	Popolare ();
	
}

void Genetic :: Quadrato ( int n, double l ) {
	
	for ( int i=0; i<n; i++) {
		double x= R->uniforme(0.,l);
		double y= R->uniforme(0.,l);
		X.push_back( {x,y} );
	}

	Popolare ();

}

void Genetic :: Circonferenza ( int n, double r ) {
	
	for ( int i=0; i<n; i++) {
		double theta= R->uniforme(0.,2.*M_PI);
		X.push_back( {r*cos(theta),r*sin(theta)} );
	}

	Popolare ();

}

void Genetic :: Print ( vector<int> Y ) {

	cout << "Stampo sul file posizioni.txt le posizioni degli oggetti allocati in un cromosoma utilizzato dal seguente algoritmo genetico" <<endl;
	
	ofstream out;
	out.open( "posizioni.txt" );
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

void Genetic :: Print ( string name, vector <int> Y ) {

	cout << "Stampo sul file " <<name << " le posizioni degli oggetti allocati in un cromosoma utilizzato dal seguente algoritmo genetico" <<endl;
	
	ofstream out;
	out.open( name );
	if ( out.is_open() ){
		for ( auto ity : Y ) {
			out << ity << "\t";
			for ( auto itx : X[ity] )
				out <<fixed <<setprecision(6) <<itx << "\t";
			out <<endl;
		}
	} else
		cerr << "PROBLEMA: Impossibile aprire " <<name << endl;
	out.close();
	
}

int Genetic :: GetGen () {

	return ngen;
	
}

int Genetic :: GetNPop () {

	return npop;
	
}

int Genetic :: GetNGeni () {

	return ngeni;
	
}

void Genetic :: PrintPopulation () {

	cout << "Stampo sul file popolazione" <<ngen << ".txt la popolazione della " <<ngen << "-esima generazione (cartella population)" <<endl; 
	
	ofstream out;
	string name= "population/popolazione" + to_string(ngen) + ".txt"; 
	out.open( name );
	if ( out.is_open() ){
		for ( int k=0; k<npop; k++ ) {
			for ( auto itc : P[k] )
				out << itc << "  ";
			out <<endl;
		}	
	} else
		cerr << "PROBLEMA: Impossibile aprire " <<name << endl;
	out.close();

}

void Genetic :: PrintSoluzione() {

	Print( Soluzione() );

}

void Genetic :: PrintSoluzione( string name ) {

	Print( name, Soluzione() );

}

double Genetic :: GetFLossGen () {

	Ordine ();
	return FLoss (P[0]);

}

double Genetic :: Get_BestHalf_FLossGen () {

	Ordine ();
	vector<double> B;
	for ( int i=0; i<npop/2; i++ )
		B.push_back( FLoss(P[i]) );
	
	return Media( B );

}

void Genetic :: Ordine () {

	sort( P.begin(), P.end(), [=]( vector<int> B, vector<int> C ) {return minFLoss(B,C);} );

}

void Genetic :: Selezione () {

	double r= R->uniforme(0.,1.), k= 5.;
	m= npop*pow(r,k);
	
	r= R->uniforme(0.,1.);
	p= npop*pow(r,k);
	//Controllo che i genitori siano diversi (m e p sono int, quindi sebbene r sia lo stesso nel cast ad int di npop*r^k si pootrebbe prendere il medesimo cromosoma come madre/padre) 
	while ( m == p ) {
		r= R->uniforme(0.,1.);
		p= npop*pow(r,k);	
	}

	F1= P[m];
	F2= P[p];
	
}

void Genetic :: Cross () {

	if ( R->uniforme(0.,1.) <= 0.7 ) {					//crossing al 70%
		int a= R->uniforme(1.,ngeni);
		for ( ; a<ngeni; a++ ) {
			for ( int j=1; j<ngeni; j++ ) {
				bool presenteM= false, presenteP= false;
				for ( int k=1; k<a; k++ ) {
					if ( P[p][j] == F1[k] )
						presenteM= true;
					if ( P[m][j] == F2[k] )
						presenteP= true;				
				}
				if ( !presenteM )
					F1[a]= P[p][j];
				if ( !presenteP )
					F2[a]= P[m][j];				
			}
		}				
	}

}

void Genetic :: Mutare () {

	//Chiamata ad una mutazione il 20% delle volte (vedi Mutare1(), Mutare2(),...)
	if ( R->uniforme (0.,1.) < 0.2 ) {
		
		int mutazione= R->uniforme(1.,5.);
		if ( mutazione == 1 )
			Mutare1 ();
		if ( mutazione == 2 )
			Mutare2 ();
		if ( mutazione == 3 )
			Mutare3 ();
		if ( mutazione == 4 )
			Mutare4 ();
	}
	
}

void  Genetic :: NewGen () {
	
	Ordine ();
	
	while ( int( P1.size() ) < npop ) {
		
		Selezione();
		Cross();
		Mutare();
		 
		P1.push_back(F1);
		P1.push_back(F2);
	
		F1.clear();
		F2.clear();
	}
	
	P.clear();
	P= P1;
	P1.clear();
	ngen++;
	
}






