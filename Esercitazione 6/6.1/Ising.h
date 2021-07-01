#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x;
double err_u,err_c,err_m,err_x;

//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk, metro;

//functions
void Input();
void Reset(int);
void Accumulate();
void Averages(int);
void Move();
void ConfFinal();
void Measure();
double PrimiVicini(int);
int Pbc(int);
double Error(double,double,int);

#endif



