#ifndef __TSP__
#define __TSP__

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include "random.h"
#include "city.h"
#include "Individual.h"
#include "Population.h"

#include "mpi.h" // Message Passing Interface(MPI) 


//Random numbers
int seed[4];
Random* rnd = new Random;

//data files
string path;

//MPI
int size, node; // node corresponds to rank (if uses rank ERROR: ambiguos use)

//parameters
int ncity,nindividual,ngen,nmigration;
double peso;

//probability
double prob_crossover;
double prob_gen_mut1,prob_gen_mut2, prob_gen_mut3, prob_gen_mut4;


//position of the first city randomly placed on a circumfernce
double xstart_city;
double ystart_city;

//cities geometry
int circle_or_square;

//starting individual
Individual starting_individual;


//starting population
Population starting_population;
Population old_pop;
Population new_pop;
Population temp;

// simulation
int index_imot, index_ifat;

//migration
int itag1, itag2, itag3, itag4, itag5, itag6;
int rk[2];


//functions
void Input(void);
void Reset(void);
void Cross_Mut(void);
void Accumulate(void);
void PrintResults(int);
void Restart(void);
void Migration(void);


//void ConfFinal(void);


#endif

