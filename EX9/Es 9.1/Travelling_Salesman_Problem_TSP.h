#ifndef __TSP__
#define __TSP__

#include <cmath>
#include "random.h"
#include "city.h"
#include "Individual.h"
#include "Population.h"


//Random numbers
int seed[4];
Random* rnd = new Random;

//parameters
int ncity,nindividual,ngen;
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


//functions
void Input(void);
void Reset(void);
void Cross_Mut(void);
void Accumulate(void);
void PrintResults(int);
void Restart(void);


//void ConfFinal(void);


#endif

