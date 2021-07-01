using namespace std;

#ifndef __Genetic__
#define __Genetic__

class Genetic {

private:
	vector<vector<double>> X;			//posizioni
	vector<vector<int>> P;				//popolazione
	vector<vector<int>> P1;				//nuova genrazione
	vector<int> A;					//vettore campione [1,2,3,...,N]
	vector<int> F1, F2;
	Random* R= new Random ();
	int m, p;
	int npop, ngeni, ngen;
	
protected:

	double Media ( vector<double> );

	double distanza ( vector<double>, vector<double> );
	double FLoss ( vector<int> );
	bool minFLoss ( vector<int> , vector<int> );				//condizione di sort di Ordine()
	
	void Popolare ();			//crea popolazione

	void Mutare1 ();			//permutazione due elemnti
	void Mutare2 ();			//movimento in avanti
	void Mutare3 ();			//permutazione di più elemeti contigui
	void Mutare4 ();			//inversione

	vector<int> Soluzione ();

public:

	Genetic ();
	~Genetic ();

	void SetGenetic ( string );
	void Quadrato ( int, double );					// ( int npart, double lato ) dove costruiamo il quadrato con un verice nell'origine 
	void Circonferenza ( int, double );				// ( int npart, double raggio ) dove consideriamo la circonferenza nell'origine
	
	void Print ( vector<int> );					// in 2D
	void Print ( string, vector<int> );				// 2D
	int GetGen ();
	int GetNPop ();
	int GetNGeni ();
	void PrintPopulation ();
	void PrintSoluzione();
	void PrintSoluzione( string );
	double GetFLossGen ();
	double Get_BestHalf_FLossGen ();
	
	void Ordine ();					//ordina la popolazione secondo bool min( vector<int>,vector<int> )
	void Selezione ();
	void Cross ();
	void Mutare ();	
	void NewGen ();
	
	
};

#endif

