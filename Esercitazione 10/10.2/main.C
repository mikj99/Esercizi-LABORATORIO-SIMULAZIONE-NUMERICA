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
#include "genetic.h"
#include "overloading.h"
#include <mpi.h>

using namespace std;

void Contaminazione ( int, int, Genetic * );

int main ( int argc, char** argv ) {
	
	//PARALLELO
	int size, rank;
	
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	//MPI_Status s1;
	//MPI_Request req;

	Genetic *G= new Genetic;
	int seed[4];
	for (int i=0; i<4; i++) {
		seed[i]= pow(i,2) + rank;
	}
	G->SetRNG( seed,2892,2587 );
	G->SetGenetic( "città.txt" );
	
	ofstream fl;
	fl.open( "FLoss_continente" + to_string( rank+1 ) + ".txt" );
	if ( rank == 0 )
		cout << "Stampo nel file FLoss_continente#.txt alcuni dati relativi alla funzione obiettivo; in paticolare:\ngenerazione | funzione obiettivo | media funzione obiettivo nella metà buona della popolazione "<<endl <<endl;

	int generazioni= 200;
	int itag0= 4, itag1= 1, itag2= 2, itag3= 3;
	int n= G->GetNGeni();
	vector<int> Scambio;
	for (int i=0; i<8; i++) 
		Scambio.push_back( 2*i + 1 );
	
	for ( int i=1; i<=generazioni; i++ ) {
		G->NewGen();
		if ( i%int(double(generazioni)/5.) == 0 ) {
			for ( auto it : Scambio ) {
				MPI_Status s0, s1,s2,s3;
				MPI_Request req0, req1, req2, req3;
				vector<int> A;
				A= G->SelectCromo (it);
				if ( rank==0 ) {
					MPI_Isend( &A[0],n,MPI_INTEGER,1,itag1,MPI_COMM_WORLD,&req0 );
					MPI_Recv( &A[0],n,MPI_INTEGER,3,itag0,MPI_COMM_WORLD,&s0 );
				}
				if ( rank==1 ) {
				 	MPI_Isend( &A[0],n,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&req1 );
					MPI_Recv( &A[0],n,MPI_INTEGER,0,itag1,MPI_COMM_WORLD,&s1 );	
				}
				if ( rank==2 ) {
				 	MPI_Isend( &A[0],n,MPI_INTEGER,3,itag3,MPI_COMM_WORLD,&req2 );
					MPI_Recv( &A[0],n,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&s2 );
				}
				if ( rank==3 ) {
				 	MPI_Isend( &A[0],n,MPI_INTEGER,0,itag0,MPI_COMM_WORLD,&req3 );
					MPI_Recv( &A[0],n,MPI_INTEGER,2,itag3,MPI_COMM_WORLD,&s3 );
				}
				G->SubCromo (it,A);	
			}
			G->Ordine();
		}
		
		fl <<fixed <<setprecision(6) <<G->GetGen() << "\t" <<G->GetFLossGen() << "\t" <<G->Get_BestHalf_FLossGen() <<endl;
		
		if ( rank==0 ) {
			if ( i%int(double(generazioni)/100.) == 0 )
				cout <<double(i)/double(generazioni)*100 << "%\t" <<flush;
			if ( i%int(double(generazioni)/10.) == 0 )
				cout <<endl;
		}
	}
	
	fl.close();
	
	vector<int> A;
	A= G->SelectCromo(0);
	MPI_Status s1,s2,s3;
	MPI_Request req1, req2, req3;
	if ( rank==1 )
	 	MPI_Isend( &A[0],n,MPI_INTEGER,0,itag1,MPI_COMM_WORLD,&req1 );
	if ( rank==2 )
		MPI_Isend( &A[0],n,MPI_INTEGER,0,itag2,MPI_COMM_WORLD,&req2 );
	if ( rank==3 ) 
		MPI_Isend( &A[0],n,MPI_INTEGER,0,itag3,MPI_COMM_WORLD,&req3 );
	if ( rank==0 ) {
		vector<int> A1=A,A2=A,A3=A;
		MPI_Recv( &A1[0],n,MPI_INTEGER,1,itag1,MPI_COMM_WORLD,&s1 );
		MPI_Recv( &A2[0],n,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&s2 );
		MPI_Recv( &A3[0],n,MPI_INTEGER,3,itag3,MPI_COMM_WORLD,&s3 );
		G->SubCromo (1,A1);
		G->SubCromo (2,A2);
		G->SubCromo (3,A3);
	
		G->Ordine();
		G->PrintSoluzione( "bestpath.txt" );
	}
	
	MPI_Finalize();
	
	return 0;
	
}


