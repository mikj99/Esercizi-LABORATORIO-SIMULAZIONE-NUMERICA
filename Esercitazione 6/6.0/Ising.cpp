#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include "Ising.h"

using namespace std;

int main() { 
	
	cout << "Come funziona il programma?" <<endl;
	cout << "1) impostare una temperatura (T) di simulazione;" <<endl;
	cout << "2) impostare una costante di campo magnetico esterno (h) agente sugli spin simulati;" <<endl;
	cout << "3) scegliere l'algoritmo di simulazione delle configurazioni;" <<endl <<endl;
	
	cout << "Temperatura? Si consiglia di mantersi circa in [0.5,3]" <<endl;
	cin >> temp;
	beta = 1.0/temp;
	
	cout << "Costante di campo magetico esterno (h)?" <<endl;
	cin >> h;
	
	cout << "Metodo con cui generare configurazioni? [M(RT)^2 -> 1; Gibbs -> 2]" <<endl;
	cin >> metro;
	if ( metro!=1 && metro!=2 ) {
		cout<< "Metodo non conforme!" <<endl;
		exit(1);
	}
		
		
	Input();
	for(int iblk=1; iblk <= nblk; ++iblk) {
		Reset(iblk);																	//Reset block averages
			for(int istep=1; istep <= nstep; ++istep) {
				Move();
				Measure();
				Accumulate();																//Update block averages
			}
		Averages(iblk);																	//Print results for current block
	cout << "Completamento: " << (iblk)*100./(double)nblk << "%" <<endl;
	}
	ConfFinal();

	return 0;
	
}


void Input() {

	cout << "Modelle 1D Classico di Ising" <<endl;
	cout << "Simulazione MonteCarlo" <<endl <<endl;
	cout << "Interazione limitata ai primi vicini" <<endl <<endl;
	cout << "Peso di Boltzmann: exp(- beta * H ), beta= 1/(k_B*T)" <<endl <<endl;
	cout << "Il programma utilizza le unità k_B=1 e mu_B=1" <<endl;

//Inizializzazione RNG															
	int seed[4];
	
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()) {
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	string nomefile= "seed.out";
	ifstream input( nomefile );
	if ( input.is_open() ){
		while ( !input.eof() ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
			}
		input.close();
	} else {
		cerr << "PROBLEM: Unable to open seed.out; i'm trying seed.in" <<endl;
		input.close();
		nomefile= "seed.in";
		input.open( nomefile );
		string property;
		if ( input.is_open() ) {
			while ( !input.eof() ) {
				input >> property;
				if ( property == "RANDOMSEED" ) {
					input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					rnd.SetRandom(seed,p1,p2);
				}
			}
			input.close();
		} else 
		cerr << "PROBLEM: Unable to open seed.in" <<endl;
	}
  
//Informazioni di base
	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> nspin;
	ReadInput >> J;
	ReadInput >> nblk;
	ReadInput >> nstep;

	cout << "Temperatura= " <<temp <<endl;
	cout << "Numero di spin coinvolti= " <<nspin <<endl;
	cout << "Inerazione di scambio= " <<J <<endl;
	cout << "Campo esterno= " <<h <<endl <<endl;
	if( metro==1 )
		cout << "Il programma utilizza il metodo di Metropolis" <<endl;
	else 
		cout << "Il programma utilizza il metodo di Gibbs" <<endl;
	cout << "Numero di blocchi= " <<nblk <<endl;
	cout << "Numbero di configurazione per blocco= " <<nstep <<endl <<endl;
	
	ReadInput.close();


//Array per le misure
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
	n_props = 4; //Number of observables

//initial configuration
	for (int i=0; i<nspin; ++i) {
    	if( rnd.uniforme(-1.,1.) > 0 ) 
    		s[i] = 1;
    	else 
    		s[i] = -1;
	}
	Measure();
	cout << "Energia iniziale= " <<walker[iu]/(double)nspin <<endl;
	
//fase equilibrazione
	for (int i=0; i<2*nstep; i++)
		Move();
	Measure();
	cout << "Energia dopo la fase di equilibrazione= " <<walker[iu]/(double)nspin <<endl;
}


void Move() {
	int o;

	for(int i=0; i<nspin; ++i) {
		o = (int)(rnd.uniforme()*nspin);											//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		if ( metro==1 ) {																//Metropolis
			double dE= 2*J*s[o]*PrimiVicini(o) + 2.*h*s[o];
			if ( dE < 0 ) 
				s[o]*= (-1);
			else
				if ( rnd.uniforme(0.,1.) < exp( (-1.)*beta*dE ) )
					s[o]*= (-1);
		} else {																			//Gibbs sampling
			double p= 1./( 1.+exp( (-2.)*beta*J*PrimiVicini(o) - 2.*h ) );
			if ( rnd.uniforme(0.,1.) < p )
				s[o]= 1.;
			else 
				s[o]= -1.;	
		}
	}
	
}

double PrimiVicini( int a ) {
	
	return ( s[Pbc(a-1)] + s[Pbc(a+1)] );

}


void Measure() {
	double u= 0.0, m= 0.0;

//cycle over spins
	for (int i=0; i<nspin; ++i) {
		u+= (-1.)*J*s[i]*s[Pbc(i+1)] - 0.5*h*( s[i]+s[Pbc(i+1)] );
		m+= s[i];
	}
	walker[iu]= u;
	walker[ic]= pow(u,2);
	walker[im]= m;
	walker[ix]= pow(m,2);

}


void Reset( int iblk ) {															//Reset block averages
   
   if( iblk == 1 ) {
		for(int i=0; i<n_props; ++i) {
           glob_av[i] = 0;
           glob_av2[i] = 0;
		}
	}

	for ( int i=0; i<n_props; ++i ) {
		blk_av[i] = 0;
	}
	blk_norm = 0;
}


void Accumulate() { 																//Update block averages

	for(int i=0; i<n_props; ++i)
		blk_av[i]+= walker[i];
	
	blk_norm+= 1.;

}


void Averages( int iblk ) {															//Print results for current block
    
	ofstream Ene, Heat, Mag, Chi;
    
	//Energy
	if ( iblk == 1 ) {
		Ene.open( "output.ene.txt" );
		Ene << "Temperatura= " <<temp <<endl;
		Ene << "Numero di spin coinvolti= " <<nspin <<endl;
		Ene << "Inerazione di scambio= " <<J <<endl;
		Ene << "Campo esterno= " <<h <<endl <<endl;
		if( metro==1 )
			cout << "Il programma utilizza il metodo di Metropolis" <<endl;
		else 
			Ene << "Il programma utilizza il metodo di Gibbs" <<endl;
		Ene << "Numero di blocchi= " <<nblk <<endl;
		Ene << "Numbero di configurazione per blocco= " <<nstep <<endl <<endl;
	} else
		Ene.open( "output.ene.txt", ios::app );
	
	stima_u = blk_av[iu]/(blk_norm*(double)nspin);
	glob_av[iu]  += stima_u;
	glob_av2[iu] += pow(stima_u,2);
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	Ene <<fixed <<setprecision(6) <<iblk << "\t" <<stima_u << "\t" <<glob_av[iu]/(double)iblk << "\t" <<err_u <<endl;
	Ene.close();

	//Heat
	if ( iblk == 1 ) {
		Heat.open( "output.heat.txt" );
		Heat << "Temperatura= " <<temp <<endl;
		Heat << "Numero di spin coinvolti= " <<nspin <<endl;
		Heat << "Inerazione di scambio= " <<J <<endl;
		Heat << "Campo esterno= " <<h <<endl <<endl;
		if( metro==1 )
			Heat << "Il programma utilizza il metodo di Metropolis" <<endl;
		else 
			Heat << "Il programma utilizza il metodo di Gibbs" <<endl;
		Heat << "Numero di blocchi= " <<nblk <<endl;
		Heat << "Numbero di configurazione per blocco= " <<nstep <<endl <<endl;
	} else 
		Heat.open( "output.heat.txt", ios::app );
	
	stima_c = pow(beta,2) * ( blk_av[ic]/(blk_norm*(double)nspin) - pow(blk_av[iu]/(blk_norm),2)/(double)nspin );
	glob_av[ic]  += stima_c;
	glob_av2[ic] += pow(stima_c,2);
	err_c= Error(glob_av[ic],glob_av2[ic],iblk);
	Heat <<fixed <<setprecision(6) <<iblk << "\t" <<stima_c << "\t" <<glob_av[ic]/(double)iblk << "\t" <<err_c <<endl;
	Heat.close();

	//Magnetization
	if ( iblk == 1 ) {
		Mag.open( "output.mag.txt" );
		Mag << "Temperatura= " <<temp <<endl;
		Mag << "Numero di spin coinvolti= " <<nspin <<endl;
		Mag << "Inerazione di scambio= " <<J <<endl;
		Mag << "Campo esterno= " <<h <<endl <<endl;
		if( metro==1 )
			Mag << "Il programma utilizza il metodo di Metropolis" <<endl;
		else 
			Mag << "Il programma utilizza il metodo di Gibbs" <<endl;
		Mag << "Numero di blocchi= " <<nblk <<endl;
		Mag << "Numbero di configurazione per blocco= " <<nstep <<endl <<endl;
	} else
		Mag.open( "output.mag.txt", ios::app );
	
	stima_m = blk_av[im]/(blk_norm*(double)nspin);
	glob_av[im]  += stima_m;
	glob_av2[im] += pow(stima_m,2);
	err_m= Error(glob_av[im],glob_av2[im],iblk);
	Mag <<fixed <<setprecision(6) <<iblk << "\t" <<stima_m << "\t" <<glob_av[im]/(double)iblk << "\t" <<err_m <<endl;
	Mag.close();

	//Magnetic susceptibility
	if ( iblk == 1 ) {
		Chi.open( "output.chi.txt" );
		Chi << "Temperatura= " <<temp <<endl;
		Chi << "Numero di spin coinvolti= " <<nspin <<endl;
		Chi << "Inerazione di scambio= " <<J <<endl;
		Chi << "Campo esterno= " <<h <<endl <<endl;
		if( metro==1 )
			Chi << "Il programma utilizza il metodo di Metropolis" <<endl;
		else 
			Chi << "Il programma utilizza il metodo di Gibbs" <<endl;
		Chi << "Numero di blocchi= " <<nblk <<endl;
		Chi << "Numbero di configurazione per blocco= " <<nstep <<endl <<endl;
	} else
		Chi.open( "output.chi.txt", ios::app );
	
	stima_x = beta*blk_av[ix]/(blk_norm*(double)nspin);
	glob_av[ix]  += stima_x;
	glob_av2[ix] += pow(stima_x,2);
	err_x= Error(glob_av[ix],glob_av2[ix],iblk);
	Chi <<fixed <<setprecision(6) <<iblk << "\t" <<stima_x << "\t" <<glob_av[ix]/(double)iblk << "\t" <<err_x <<endl;
	Chi.close();

}


void ConfFinal() {

	cout << "Stampo la configurazione finale nel file config.final.txt" <<endl <<endl;
	
	ofstream WriteConf;
	WriteConf.open("config.final.txt");
	for (int i=0; i<nspin; ++i) {
		WriteConf << s[i] << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}

int Pbc( int i ) { 																	//Algorithm for periodic boundary conditions

	if( i >= nspin)  
	i+= (-1.)*nspin;
	else 
		if(i < 0) 
			i+= nspin;
    
    return i;

}

double Error( double sum, double sum2, int iblk ) {
	
	if( iblk==1 ) 
		return 0.0;
	else 
		return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2) )/(double)(iblk-1));

}




