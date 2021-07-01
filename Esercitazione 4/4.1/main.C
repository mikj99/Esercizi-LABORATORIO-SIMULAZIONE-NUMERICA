#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

//parametri ed osservabili
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

//medie
double acc,att;

//connfigurazione
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

//stati termodinamici
int npart;
double energy,temp,vol,rho,box,rcut;

//simulazione
int nstep, iprint, seed;
double delta;
bool first;

//prototipi funzioni
void Input();
void LoadConfig();
void CreateConfig();
void Rescale();
void Partenza();
void Termalizzare();
void Move();
void ConfFinal();
void ConfXYZ(int);
void Measure(int);
double Force(int, int);
double Pbc(double);

//main
int main() { 
  
	cout << "Come funziona il programma:" <<endl;
	cout << "1) essere in possesso di una buona configurazione precedente di partenza (nel caso crearla);" <<endl;
	cout << "2) eseguire la simulazione (nel caso rieseguendo il programma);" <<endl <<endl;
	cout << "Preparo una configurazione di partenza adeguata exnovo? [si 1/ no 0]" <<endl;
	cin >> first;
  	cout <<endl;
  
	Input();
	Partenza();
	if ( first ) {
		CreateConfig();
	} else {
		LoadConfig();
	}
	Termalizzare();
	
	cout << "Simulazione:" <<endl;
	int nconf = 1;
	for ( int istep=1; istep <= nstep; ++istep ) {
		Move();
		if ( istep%iprint == 0 ) 
			cout << "Numero di time-steps: " << istep << endl;
		if( istep%10 == 0 ) {
			Measure(nconf);
			//ConfXYZ(nconf);
			nconf++;
		}
	
	}
	ConfFinal();
	
	return 0;
	
}

// Funzioni

void Input() {
  
	ifstream ReadInput;
	double ep, ek, pr, et, vir;

	seed = 1;
	srand(seed);
  
	ReadInput.open("input.dat");

	ReadInput >> temp;
	ReadInput >> npart;
	ReadInput >> rho;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	ReadInput >> rcut;
	ReadInput >> delta;
	if ( first ) {
		nstep= 1500;
	} else {
		ReadInput >> nstep;
	}
	ReadInput >> iprint;

 	ReadInput.close();

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	n_props = 4; //Number of observables
   	
	cout << "Fluido classico di Lennard-Jones" << endl;
	cout << "Simulazione di dinamica molecolare in insieme NVE" << endl << endl;
	cout << "Potenziale di interazione v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Il programma usa unità di Lennard-Jones" << endl;
   	cout << "Grandezza del box di simulazione = " << box << endl;
   	cout << "Numero di particelle = " << npart << endl;
   	cout << "Densità delle particelle = " << rho << endl;
   	cout << "Volume del box di simulazione = " << vol << endl;
	cout << "Il programma integra le equazioni di Newton von l'algoritmo di Verlet" << endl;
	cout << "Time step = " << delta << endl;
	cout << "Numero di steps = " << nstep << endl << endl;
   	
	return;
}

void LoadConfig() {
	ifstream ReadConf;
  	ReadConf.open ("old.0");
  	if ( !ReadConf.is_open() ) {
  		cerr << "Configurazione precedente non esistente; creane una!" <<endl;
  		exit(1);
 	} else {	
  		for (int i=0; i<npart; ++i) {
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] *= box;   	 			
			yold[i] *= box;
			zold[i] *= box;
 		}
 	}
	ReadConf.close();
		
	ReadConf.open ("old.final");
	if ( !ReadConf.is_open() ) {
 		cerr << "Configurazione di partenza collegata ad una precedente non esistente; creane una!" <<endl;	
  		exit(1);
  	} else {
		for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i];
    			x[i] *= box;
			y[i] *= box;
			z[i] *= box;
  		}
  	}
  	ReadConf.close();
	
	//velocità:
	for (int i=0; i<npart; ++i) {
		vx[i] = Pbc ( (x[i]-xold[i])/delta );
		vy[i] = Pbc ( (y[i]-xold[i])/delta );
		vz[i] = Pbc ( (y[i]-xold[i])/delta );
	}
	Rescale();

}

void CreateConfig() {
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
	if ( first ) 
		return;
	
	int a= 1000;	
	for ( int k=1; k<11; k++ ) {
		cout << "Ripartenza " <<k << ", " <<a << " steps." <<endl; 
		LoadConfig();
		for ( int i=0; i<a; i++ )
			Move();
		ConfFinal();
	}	

}

void Termalizzare() {
	if ( first ) 
		return;
	
	cout<< "Termalizzazione: " <<nstep/2 << " steps" <<endl;
	for ( int i=0; i< nstep/2; i++ )
		Move();

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

void Measure( int nconf ){
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream out;

	if ( nconf==1 ) {
		out.open("out.txt");
		cout << "Stampa delle variabili termodinamiche di interesse delle simulazione in out.txt secondo K|P|Etot|T" <<endl;
	} 
	else
		out.open("out.txt",ios::app);

	v = 0.0; 																		//reset observables
	t = 0.0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc( xold[i] - xold[j] );
			dy = Pbc( yold[i] - yold[j] );
			dz = Pbc( zold[i] - zold[j] );

			dr = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) );

			if ( dr<rcut ) {
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);								//Potenzial energy
				v += vij;
			}
		}          
	}

//Kinetic energy
	for (int i=0; i<npart; ++i) 
		t += 0.5 * (pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
   
	stima_pot = v/(double)npart;
	stima_kin = t/(double)npart;
	stima_temp = (2.0 / 3.0) * t/(double)npart;
	stima_etot = (t+v)/(double)npart;

	out <<fixed <<setprecision(6) <<stima_kin << "\t" <<stima_pot << "\t" << stima_etot << "\t" << stima_temp << endl;
	out.close();
	
	return;
}


void ConfFinal() {
  
	ofstream WriteConf;
		
	cout << "Stampa della configurazione finale in olf.final " << endl;
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

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc( double r ){												//Periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


