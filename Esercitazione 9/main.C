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

using namespace std;

int main () {

	Genetic* G= new Genetic ();

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
	
	cout << "Scegliere per quante generazioni (numero intero) eseguire l'evoluzione" <<endl;
	int generazioni;
	cin >>generazioni;
	cout <<endl;
	
	if ( spazio=="quadrato" )
		G->Quadrato( n, 10. );
	if ( spazio=="circonferenza" )
		G->Circonferenza( n, 1. );
	if ( spazio=="file esterno" )
		G->SetGenetic( nome );

	//G->PrintPopulation();
	
	ofstream fl;
	fl.open( "FLoss_" + spazio + ".txt" );
	cout << "Stampo nel file FLoss_" <<spazio << ".txt alcuni dati relativi alla funzione obiettivo; in paticolare:\ngenerazione | funzione obiettivo | media funzione obiettivo nella metà buona della popolazione "<<endl <<endl;
	
	cout << "Eseguo per " <<generazioni << " generazioni:" <<endl;
	for ( int i=1; i<=generazioni; i++ ) {
		
		G->NewGen();
		//G->PrintPopulation();
		
		fl <<G->GetGen() << "\t" <<G->GetFLossGen() << "\t" <<G->Get_BestHalf_FLossGen() <<endl;
		
		if ( i%int(double(generazioni)/100.) == 0 )
			cout <<double(i)/double(generazioni)*100 << "%\t" <<flush;
		if ( i%int(double(generazioni)/10.) == 0 )
			cout <<endl;

	}
	cout<<endl;
	
	fl.close();
	
	G->PrintSoluzione( "bestpath_" + spazio + ".txt" );

	delete G;					//per far attivare volontariamente distruttore di G
	   
	return 0;
}



