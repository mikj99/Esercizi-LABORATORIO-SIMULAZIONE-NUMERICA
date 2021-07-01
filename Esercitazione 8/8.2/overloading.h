#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iomanip>

using namespace std;

//overloading somma +

template <typename T> vector <T> operator+ ( const vector <T> &a, const vector <T> &b ) {

	assert ( a.size() == b.size() );
	
	vector <T> tot ( a.size() );
	for (unsigned int i=0; i<a.size(); i++ )
		tot[i]= a[i]+b[i];

	return tot;

}

//overloading differenza -

template <typename T> vector <T> operator- ( const vector <T> &a, const vector <T> &b ) {

	assert ( a.size() == b.size() );
	
	vector <T> tot ( a.size() );
	for (unsigned int i=0; i<a.size(); i++ )
		tot[i]= a[i]-b[i];

	return tot;

}

//overloading prodotto scalare *

template <typename T> T operator* ( const vector <T> &a, const vector <T> &b ) {

	assert ( a.size() == b.size() );
	
	T tot {0};
	for (unsigned int i=0; i<a.size(); i++ )
		tot+= a[i]*b[i];

	return tot;

}

//overloading prodotto per scalare * (scalare,vettore)

template <typename T> vector <T> operator* ( T a, const vector <T> &b ) {
	
	vector <T> tot ( b.size() );
	for (unsigned int i=0; i<b.size(); i++ )
		tot[i]= a*b[i];

	return tot;

}

//overloading prodotto per scalare * (vettore,scalare)

template <typename T> vector <T> operator* ( const vector<T> &a, T b ) {
	
	vector <T> tot ( a.size() );
	for (unsigned int i=0; i<a.size(); i++ )
		tot[i]= a[i]*b;

	return tot;

}

//overloading divisione per scalare /

template <typename T> vector <T> operator/ ( const vector<T> &a, T b ) {
	
	vector <T> tot ( a.size() );
	for (unsigned int i=0; i<a.size(); i++ )
		tot[i]= a[i]/(double)b;

	return tot;

}

//overloading somma consecutiva +=

template <typename T> void operator+= ( vector <T> &a, const vector <T> &b ) {

	assert ( a.size() == b.size() );
	
	for (unsigned int i=0; i<a.size(); i++ )
		a[i]+= b[i];

}

//overloading somma consecutiva -=

template <typename T> void operator-= ( vector <T> &a, const vector <T> &b ) {

	assert ( a.size() == b.size() );
	
	for (unsigned int i=0; i<a.size(); i++ )
		a[i]-= b[i];

}

//stampa 

template <typename T> void Stampa ( const vector <T> &a ) {
	
	//cout << "Stampa vector" <<endl;
    for ( auto it:a ) 
		cout <<it << "  " ;
    cout <<endl;
	//cout << "Termine stampa vector" <<endl;

}

//stampa con setprecision

template <typename T> void Stampa ( const vector <T> &a, int b ) {
	
    for ( auto it:a ) 
		cout <<setprecision(b) <<it << "  " ;
    cout <<endl;

}

//stampa se l'operatore vuole stampare

template <typename T> void Stampa ( const vector <T> &a, bool b ) {
	
    if ( b ) {
		for ( auto it:a ) 
			cout <<setprecision(b) <<it << "  " ;
    	cout <<endl;
	}

}

//stampa file+setprecision

template <typename T> void Stampa ( const vector <T> &a, int b, ofstream &out ) {
	
    for ( auto it:a ) 
		out <<fixed <<setprecision(b) <<it << "\t" ;
    out <<endl;

}







