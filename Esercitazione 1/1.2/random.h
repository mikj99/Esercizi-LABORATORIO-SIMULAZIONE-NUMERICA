#ifndef __Random__
#define __Random__

class Random {

private:
	int m1;
	int m2;
	int m3;
	int m4;
	int l1;
	int l2;
	int l3;
	int l4;
	int n1;
	int n2;
	int n3;
	int n4;

protected:

public:

	Random();
	~Random();

	void SetRandom( int *, int, int );
	void SaveSeed();
	double uniforme();
	double uniforme( double, double );					// ( double min, double max )
	double Gauss( double, double );					// ( double mean, double sigma )
	double exp ( double );					// ( double lambda )
	double Lorentz ( double, double );					// ( double Gamma, double mean )

};

#endif

