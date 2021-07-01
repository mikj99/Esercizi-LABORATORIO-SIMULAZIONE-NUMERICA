#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <string>
#include <vector>
#include "Monte_Carlo_NVT.h"

using namespace std;

//stato simulazione:
string state;


int main() { 
	
	cout << "\nPrima della simulazione scegliere quale fase simulare [solid/liquid/gas]:" <<endl;
	cin >> state;
	if ( state!="solid" && state!="liquid" && state!="gas" ) {
		cout<< "Fase non conforme!" <<endl;
		exit(1);
	}
	cout <<endl;
	
	Input();
//	int nconf = 1;
	for(int iblk=1; iblk <= nblk; ++iblk) {
		Reset(iblk);
		for(int istep=1; istep <= nstep; ++istep) {
			Move();
			Measure();
			Accumulate();
//			if( istep%10 == 0 ) {
//				ConfXYZ(nconf);
//				nconf += 1;
//			}
		}
		Averages(iblk);
	}
//	ConfFinal();

	return 0;
}


void Input()	{
	
	ifstream ReadInput,ReadConf;

	cout << "Fluido Classico di Lennard-Jones" << endl;
	cout << "Simulazione Monte Carlo" << endl << endl;
	cout << "Potenziale di interazione: v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Peso di Boltzmann exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
	cout << "Il programma utilizza le unità di Lennard-Jones" << endl;

//Read seed for random numbers
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

  
//Read input informations
	ReadInput.open( "input." + state );
	
	ReadInput >>temp;
	beta= 1.0/temp;
	ReadInput >>npart;
	ReadInput >>rho;
	vol= (double)npart/rho;
	box= pow(vol,1./3.);
	ReadInput >>rcut;	
	
	//Tail corrections for potential energy and pressure
  	vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  	ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  	
	ReadInput >>delta;
	ReadInput >>nblk;
	ReadInput >>nstep;	
  	
  	
  	cout << "Temperatura= " <<temp <<endl;
	cout << "Numero di particelle= " <<npart <<endl;
	cout << "Densità di particelle= " <<rho <<endl;
	cout << "Volume del box di simulazione= " <<vol <<endl;
	cout << "Grandezza del box di simulazione= " <<box <<endl;
	cout << "Raggio di cutoff del potenziale di interazione= " <<rcut <<endl <<endl;
  	cout << "Correzione di coda dell'enegia di potenziale= " <<vtail << endl;
  	cout << "Correzione di coda del viriale = " << ptail << endl; 
	cout << "Il programma utilizza un movimento secondo l'algoritmo di Metropolis con una proposta di transizione uniforme" << endl;
	cout << "Parametro di movimento= " << delta << endl;
	cout << "Numero di blocchi= " << nblk << endl;
	cout << "Numero di step per blocco= " << nstep << endl << endl;
	ReadInput.close();


//Prepare arrays for measurements
	iv = 0;
	iw = 1;
	n_props = 2;
	
//measurement of g(r)
	igofr = 2;
	nbins = 100;
	n_props+= nbins;
	bin_size = (box/2.)/(double)nbins;

//Read initial configuration
	cout << "Lettura della configurazione iniziale secondo il file config.0" << endl << endl;
 	ReadConf.open("config.0");
	for ( int i=0; i<npart; ++i ) {
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = Pbc( x[i] * box );
		y[i] = Pbc( y[i] * box );
		z[i] = Pbc( z[i] * box );
	}
	ReadConf.close();
  
//Equilibrazione + condizione iniziale
	for (int i=0; i<int(2.*nstep*nblk/100.); i++)
		Move();
	Measure();
	cout << "Sistama equilibrato per " <<int(2.*nstep*nblk/100.) << " simulazioni;" <<endl; 
	cout << "\t- Energia potenziale iniziale (con correzione di coda) = " << walker[iv]/(double)npart + vtail << endl;
	cout << "\t- Viriale                     (con correzione di coda) = " << walker[iw]/(double)npart + ptail << endl;
	cout << "\t- Pressione                   (con correzione di coda) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;

}


void Move() {
	
	int o;
	double energy_old, energy_new;
	double xold, yold, zold, xnew, ynew, znew;


	for ( int i=0; i<npart; ++i ) {
		o = (int)(rnd.uniforme()*npart);			 //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)

		//Old
		xold = x[o];
		yold = y[o];
		zold = z[o];
		energy_old = Boltzmann(xold,yold,zold,o);

		//New
		xnew = Pbc( x[o] + rnd.uniforme( (-1.)*delta*0.5, delta*0.5) );
		ynew = Pbc( y[o] + rnd.uniforme( (-1.)*delta*0.5, delta*0.5) );
		znew = Pbc( z[o] + rnd.uniforme( (-1.)*delta*0.5, delta*0.5) );
		energy_new = Boltzmann(xnew,ynew,znew,o);

		//Metropolis test
		if ( rnd.uniforme(0.,1.) < exp( (-1.)*beta*(energy_new-energy_old) ) )  {
			x[o] = xnew;
			y[o] = ynew;
			z[o] = znew;
			accepted++;
		}	
		attempted++;
	}

}

double Boltzmann( double xx, double yy, double zz, int ip ) {

	double ene= 0.0;
	double dr;

	for ( int i=0; i<npart; ++i ) {
    		if( i != ip ) {
			dr = sqrt( pow(Pbc(xx - x[i]),2) + pow(Pbc(yy - y[i]),2) + pow(Pbc(zz - z[i]),2) );
			if( dr < rcut ) {
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}

	return 4.*ene;
}

void Measure() {
	
	double v = 0.0, w = 0.0;
	double dr;

	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) 
		walker[k]=0.0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i) {
		for (int j=i+1; j<npart; ++j) {
     			dr = sqrt( pow(Pbc(x[i] - x[j]),2) + pow(Pbc(y[i] - y[j]),2) + pow(Pbc(z[i] - z[j]),2) );			// distance i-j in pbc
			walker[int(dr/bin_size)+igofr]+= 2;
	
			if ( dr<rcut ) {
				// contribution to energy and virial
				v += 1.0/pow(dr,12) - 1.0/pow(dr,6);
				w += 1.0/pow(dr,12) - 0.5/pow(dr,6);
			}
		}          
	}

	walker[iv] = 4.*v;
	walker[iw] = w*48./3.;

}


void Reset( int iblk ) {
   
	if(iblk == 1) {
		for(int i=0; i<n_props; ++i) {
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;

}


void Accumulate() { //Update block averages

	for(int i=0; i<n_props; ++i)
		blk_av[i] += walker[i];
		
	blk_norm++;

}


void Averages( int iblk ) {							//Print results for current block
   
	double dv;
	ofstream out, Gout;
    
	cout << "Numero di blocco " << iblk << endl;
	cout << "Rateo di accetazione= " << accepted/attempted << endl << endl;
    
  	if ( iblk == 1 ) {
		out.open( "out." + state + ".txt" ); 
	} else {
		out.open( "out." + state + ".txt", ios::app );
	}

    
	stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail;					//per particella
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
	stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol;
	glob_av[iw] += stima_pres;
	glob_av2[iw] += stima_pres*stima_pres;
	err_press=Error(glob_av[iw],glob_av2[iw],iblk);

	for (int i=igofr; i<igofr+nbins; i++) {
		dv= 4.*M_PI*(pow(bin_size*i,3)-pow(bin_size*(i-1),3))/3.;
		blk_av[i]/= rho*(double)npart*dv*(double)nstep;
		glob_av[i]+= blk_av[i];
		glob_av2[i]+= pow(blk_av[i],2);
	}

    	out <<fixed <<setprecision(6) <<iblk << "\t" <<stima_pot << "\t" <<glob_av[iv]/(double)iblk << "\t" <<err_pot << "\t" << stima_pres << "\t" << glob_av[iw]/(double)iblk << "\t" <<err_press <<endl;
	if ( iblk == nblk ) {
		Gout.open ( "out.g(x)." + state + ".txt" ); 	
		for (int i=igofr; i<igofr+nbins; i++) {
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			Gout <<fixed <<setprecision(6) <<bin_size*(i-1) << "\t" <<bin_size*(i) << "\t" <<glob_av[i]/(double)iblk << "\t" <<err_gdir <<endl;
		}
		Gout.close();
	}
	cout << "----------------------------" << endl << endl;

	out.close();
	
}


void ConfFinal() {
	
	ofstream WriteConf;
	cout << "Stampa della configurazione finale nel file config.final " << endl << endl;
	WriteConf.open("config.final");
	for ( int i=0; i<npart; ++i ) {
		WriteConf <<x[i]/box << "   " <<y[i]/box << "   " <<z[i]/box << endl;
	}
	
	WriteConf.close();

	rnd.SaveSeed();
	
}

void ConfXYZ( int nconf ) {			//Write configuration in .xyz format
  
  	ofstream WriteXYZ;
	WriteXYZ.open("frames/config_" + to_string(nconf) + state +".xyz");
	WriteXYZ << npart << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	
	WriteXYZ.close();

}

double Pbc( double r ) {		//Algorithm for periodic boundary conditions with side L=box

    return r - box * rint(r/box);
    
}

double Error(double sum, double sum2, int iblk) {

	if( iblk == 1 ) 
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
		
}

