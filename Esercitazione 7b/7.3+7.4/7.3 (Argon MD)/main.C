#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <string>

using namespace std;

//parametri ed osservabili
const int m_props=1000;
int n_props, iv,ik,it,ie,iw,igofr;
double bin_size, nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot, stima_kin, stima_temp, stima_etot, stima_pres, err_pot, err_kin, err_temp, err_etot, err_press, err_gdir;

//connfigurazione
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

//stati termodinamici
int npart;
double energy,temp,vol,rho,box,rcut;
double ptail,vtail;

//simulazione
int nstep, nmeasure, nblk, seed;
double delta;
string state;

//prototipi funzioni
void Input();
void CreateConfig();
void Rescale();
void Termalizzare();
void Partenza();
void Reset( int );
void Accumulate();
void Averages( int );
void Move();
void ConfFinal();
void ConfXYZ(int);
void Measure();
double Force(int, int);
double Pbc(double);
double Error( double, double, int );


//main
int main() { 
  
	cout << "\nPrima della simulazione scegliere quale fase simulare [solid/liquid/gas]:" <<endl;
	cin >> state;
	if ( state!="solid" && state!="liquid" && state!="gas" ) {
		cout<< "Fase non conforme!" <<endl;
		exit(1);
	}
	cout <<endl;
	
	Input();
	CreateConfig();
	Partenza();
	Termalizzare();

	cout << "Simulazione:" <<endl;
	for(int iblk=1; iblk <= nblk; ++iblk) {
		Reset(iblk);
		for(int istep=1; istep <= nstep; ++istep) {
			Move();
			if( istep%nmeasure == 0 ) {
				Measure();
				Accumulate();
			}

		}
		Averages(iblk);
	}
	
	return 0;
	
}

// Funzioni

void Input() {

	seed = 1;
	srand(seed);
  
  	ifstream ReadInput;
	ReadInput.open( state+".dat" );

	ReadInput >> temp;
	ReadInput >> npart;
	ReadInput >> rho;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	ReadInput >> rcut;
	
	//Tail corrections for potential energy and pressure
  	vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  	ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
	
	ReadInput >> delta;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput >> nmeasure;

 	ReadInput.close();

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	iw = 4; //virial
	n_props = 5; //Number of observables
	
	//measurement of g(r)
	igofr = 5;
	nbins = 100;
	n_props+= nbins;
	bin_size = (box/2.)/(double)nbins;
   	
	cout << "Fluido classico di Lennard-Jones" << endl;
	cout << "Simulazione di dinamica molecolare in insieme NVE" << endl << endl;
	cout << "Potenziale di interazione v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Il programma usa unità di Lennard-Jones" << endl;
   	cout << "Grandezza del box di simulazione = " << box << endl;
   	cout << "Numero di particelle = " << npart << endl;
   	cout << "Densità delle particelle = " << rho << endl;
   	cout << "Volume del box di simulazione = " << vol << endl;
   	cout << "Raggio di cutoff del potenziale di interazione= " <<rcut <<endl <<endl;
  	cout << "Correzione di coda dell'enegia di potenziale= " <<vtail << endl;
  	cout << "Correzione di coda del viriale = " << ptail << endl; 
	cout << "Il programma integra le equazioni di Newton con l'algoritmo di Verlet" << endl;
	cout << "Time step = " << delta << endl;
	cout << "Numero di blocchi= " << nblk << endl;
	cout << "Numero di step per blocco= " << nstep/nmeasure << endl << endl;
   	
	return;
}

void CreateConfig() {
		cout << "Creo configurazione di partenza." <<endl <<endl;
		ifstream ReadConf;
		ReadConf.open("config.0");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] *= box;
			y[i] *= box;
    			z[i] *= box;
    		}
 		ReadConf.close();

		//velocità:
		for (int i=0; i<npart; ++i) {
			vx[i] = rand()/double(RAND_MAX) - 0.5;
			vy[i] = rand()/double(RAND_MAX) - 0.5;
			vz[i] = rand()/double(RAND_MAX) - 0.5;
		}
		Rescale();

}

void Rescale() {		
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i) {
			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) 
			sumv[idim] /= (double)npart;
 		
 		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] -= sumv[0];
			vy[i] -= sumv[1];
			vz[i] -= sumv[2];
			sumv2 += pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta); 
   		}
}

void Partenza() {

	for ( int k=1; k<11; k++ ) {
		cout << "Ripartenza " <<k << ", 1000 steps." <<endl; 
		for ( int i=0; i <1000; ++i )
			Move();
		Rescale();
	}	
	cout <<endl;

}

void Termalizzare() {
	
	cout<< "Termalizzazione: " <<nstep/2 << " steps." <<endl;
	for ( int i=1; i<=nstep/2; i++ ) {
		Move();
		if ( i%(nstep/20) == 0 )
			cout <<i/(nstep/200) << "%\t" <<flush;
	}
	cout <<endl;
	
	Measure();
	cout << "Condizioni iniali della simulazione:" <<endl;
	cout << "\t- Energia potenziale iniziale (con correzione di coda) = " << walker[iv]/(double)npart + vtail << endl;	
	cout << "\t- Energia cinetica iniziale                            = " << walker[ik]/(double)npart << endl;
	cout << "\t- Energia totale inizle                                = " << walker[ie]/(double)npart << endl;
	cout << "\t- Temperatura inizle                                   = " << walker[it]/(double)npart << endl;
	cout << "\t- Pressione iniziale          (con correzione di coda) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;

}



void Move() {
	
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for (int i=0; i<npart; ++i) {
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
  }

	for(int i=0; i<npart; ++i) {													//Algoritmo di Verlet

	xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );						//Posizioni di Verlet
    	ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    	znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    	vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);									//Velocità di Verlet
    	vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    	vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    	xold[i] = x[i];
    	yold[i] = y[i];
    	zold[i] = z[i];
    	x[i] = xnew;
    	y[i] = ynew;
    	z[i] = znew;
	}
	
	return;
}

double Force( int ip, int idir ) {												//Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

  for ( int i=0; i<npart; ++i ) {
    if( i!=ip ) {
      dvec[0] = Pbc( x[ip] - x[i] );
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = pow(dvec[0],2) + pow(dvec[1],2) + pow(dvec[2],2);
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  
  return f;
}

void Measure() {
	int bin;
	double v=0.0, w=0.0, t=0.0, dr;
	
	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) 
		walker[k]=0.0;


	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j) {
			dr = sqrt( pow(Pbc( xold[i] - xold[j] ),2) + pow(Pbc( yold[i] - yold[j] ),2) + pow(Pbc( zold[i] - zold[j] ),2) );
			
			//Potenzial energy
			if ( dr<rcut ) {
				v+= 4.0/pow(dr,12) - 4.0/pow(dr,6);						//Potenzial energy
				w+= 1.0/pow(dr,12) - 0.5/pow(dr,6);
			}
			
			//G(x)
			walker[int(dr/bin_size)+igofr]+= 2;
			
		}          
	}

	//Kinetic energy
	for (int i=0; i<npart; ++i) 
		t += 0.5 * (pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
   
	walker[iv] = v;
	walker[ik] = t;
	walker[it] = (2.0 / 3.0) * t;
	walker[ie] = (t+v);
	walker[iw] = w*48./3.;
	
	return;
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

}


void Accumulate() { //Update block averages

	for(int i=0; i<n_props; ++i)
		blk_av[i] += walker[i];
		
	blk_norm++;

}

void Averages( int iblk ) {							//Print results for current block
   
	double dv;
	ofstream out, Gout;
    
	cout << "Esecuzione blocco " << iblk << endl;
    
  	if ( iblk == 1 ) {
		out.open( "out." + state + ".txt" ); 
	} else {
		out.open( "out." + state + ".txt", ios::app );
	}

    
	stima_pot = blk_av[iv]/blk_norm/(double)npart +vtail;
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
	
	stima_kin = blk_av[ik]/blk_norm/(double)npart;
	glob_av[ik] += stima_kin;
	glob_av2[ik] += pow(stima_kin,2);
	err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
	
	stima_temp = blk_av[it]/blk_norm/(double)npart;
	glob_av[it] += stima_temp;
	glob_av2[it] += pow(stima_temp,2);
	err_temp=Error(glob_av[it],glob_av2[it],iblk);
	
	stima_etot = blk_av[ie]/blk_norm/(double)npart;
	glob_av[ie] += stima_etot;
	glob_av2[ie] += pow(stima_etot,2);
	err_temp=Error(glob_av[ie],glob_av2[ie],iblk);
	
	stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol;
	glob_av[iw] += stima_pres;
	glob_av2[iw] += stima_pres*stima_pres;
	err_press=Error(glob_av[iw],glob_av2[iw],iblk);

	for (int i=igofr; i<igofr+nbins; i++) {
		dv= 4.*M_PI*(pow(bin_size*(i-3),3)-pow(bin_size*(i-4),3))/3.;
		blk_av[i]/= rho*(double)npart*dv*(double)nstep/nmeasure;
		glob_av[i]+= blk_av[i];
		glob_av2[i]+= pow(blk_av[i],2);
	}

    	out <<fixed <<setprecision(6) <<iblk << "\t" <<stima_pot << "\t" <<glob_av[iv]/(double)iblk << "\t" <<err_pot << "\t"
    						     <<stima_kin << "\t" << glob_av[ik]/(double)iblk << "\t" <<err_kin << "\t"
    						     <<stima_temp << "\t" << glob_av[it]/(double)iblk << "\t" <<err_temp << "\t"
    						     <<stima_etot << "\t" << glob_av[ie]/(double)iblk << "\t" <<err_etot << "\t"
    						     <<stima_pres << "\t" << glob_av[iw]/(double)iblk << "\t" <<err_press <<endl;
    						     
	if ( iblk == nblk ) {
		Gout.open ( "out.g(x)." + state + ".txt" ); 	
		for (int i=igofr; i<igofr+nbins; i++) {
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			Gout <<fixed <<setprecision(6) <<bin_size*(i-4) << "\t" <<bin_size*(i-3) << "\t" <<glob_av[i]/(double)iblk << "\t" <<err_gdir <<endl;
		}
		Gout.close();
	}

	out.close();
	
}

void ConfFinal() {
  
	ofstream WriteConf;
		
	cout << "Stampa della configurazione finale in old.final " << endl;
	WriteConf.open("old.final");
	for (int i=0; i<npart; ++i){		
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}
	WriteConf.close();
    	
	cout << "Stampa della configurazione precedente in old.0 " << endl << endl;
	WriteConf.open("old.0");
	for (int i=0; i<npart; ++i){
		WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
	}
	WriteConf.close();
  
	return;
}

void ConfXYZ( int nconf ) {												//Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + state + ".xyz");
  WriteXYZ << npart << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc( double r ){												//Periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk) {

	if( iblk == 1 ) 
		return 0.0;
	else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
		
}



