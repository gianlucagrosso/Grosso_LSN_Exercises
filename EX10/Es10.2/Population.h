#ifndef __Population_h__
#define __Population_h__

#include "Individual.h"


class Population {
	public:	
		void NewIndividual (Individual);

		Individual GetIndividual(int);
		long unsigned int GetNIndividual(); //return the number of individuals present in the population
		void Setcity(int,int,double,double,int);

		bool Check(); //check function which verify if every individuals in the population fulfils the bonds

		void Clear() {m_pop.clear();}; //cancel the population

		void StartingPop(); //initialize the starting population
	
		int Selection(int,Random*,double); //selection operator

		void SortPop(); // sort the population as a function of the increasing Cost function of each individuals 
	
		void PrintPop(const char* Filename); 
		void PrintBest(int,const char* Filename,int); //Print the best individuals in the ordered population
		void Print_The_Best_Path(int,const char* Filename); //Print the best individuals in absolute
		void PrintHalfAverage(const char* Filename,int,int); //print the averaged length value on the best half of the population

		void Mutation(int,int,int,Random*); //Genetic mutation operators
		void CrossOver(int,int,int); //Crossover operator	
		
	private:
		vector <Individual> m_pop; // population: Individuals' contenitor
};
#endif
