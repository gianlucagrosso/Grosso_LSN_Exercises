#ifndef __VMC__
#define __VMC__


#include "random.h"

//Random numbers
int seed[4];
Random rnd;

//variational parameters
double mu,sigma;
double temp , beta , energy = 0. , energy_err = 0.;

//global parameters
int iQM1D , iPrint;

//position
double x_pos;

//parameters, observables
const int m_props=5000;
int n_props, ih, iw; //hamiltonian and histogram indexes
double bin_size,nbins, histogram_start, histogram_end;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_ham,err_ham, stima_PsiT2 , err_PsiT2;
std::ofstream Ham, PsiT_2;

// VMC simulation 
int nstep, nblk , Nblk_appo;
double deltaVMC , VMCaccepted, VMCattempted;;
double start_position;

//SA simulation
double deltaSA , SAaccepted , SAattempted;
int NSA_step;
double en_min , mu_min , sigma_min, en_min_err;
std::vector<double> delta={2.8, 1.6, 1.5, 1.0, 0.7, 0.4, 0.235, 0.125, 0.06, 0.01}; //set of temprature and delta to have SA acceptance 0.5
std::vector<double> T = {50, 18, 7, 2.5, 0.9 , 0.3 , 0.12 , 0.05 , 0.015 , 0.006};


//functions


void Input(void);

//Variational Monte Carlo
void VMCReset(int);
void VMCAccumulate(void);
void VMCAverages(int);
void VMCMove(void);
void VMCMeasure(void);
void VMCConfFinal(void);

//Simulated Annealing
void SAAverages(int);
void SAMove(void);
void SAReset(void);
void Energy(int);
void FindDelta(void);


double Error(double,double,int);
double Psi(double, int);
double DerSecPsi(double , int); 
double Pot(double , int);
double Pbc(double);

#endif

