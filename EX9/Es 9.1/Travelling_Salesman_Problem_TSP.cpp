
#include "Travelling_Salesman_Problem_TSP.h"
#include "funzioni.h"


using namespace std;

int main (int argc, char *argv[]){

   Input(); //initialitation
   
   for(int igen=1; igen<=ngen; igen++)
   {
	for(int jInd=0; jInd<nindividual-1; jInd+=2)
	{
		Reset(); //Reset and riinitialize temp population
		Cross_Mut(); //Do the crossover and the mutation
		Accumulate(); //Accumulate new individuals in new population	
	}
	new_pop.SortPop(); //Sort the new population
	if (!starting_population.Check()) cerr << "PROBLEM: The population doesn't fulfil the bonds of the problem" << endl << endl; //check

	PrintResults(igen); 	
	Restart(); //Rearrange the populations for the next generational cycle	
   }
			  
   rnd->SaveSeed();
   return 0;
}


void Input(void){

  ifstream ReadInput;

  //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   cout << "This program solve the Travelling salesman Problem with a Genetical Algorithm" << endl;
   cout << "The bonds are: the salesman must visit one and only one every city and must be back to the fisrt city at the end of the path!" << endl << endl;


  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> ncity;
  cout << "The path is made of " << ncity << " cities" << endl;

  ReadInput >> nindividual;
  cout << "Each poplation is made of " << nindividual << " individuals" << endl <<endl;

  ReadInput >> ngen;
  cout << "The generations (new populations) produced by the algorithm are : " << ngen << endl;

  ReadInput >> xstart_city; 
  ReadInput >> ystart_city;

  ReadInput >> circle_or_square;

//The starting and returning city of each paths
  if (circle_or_square==1)  //city randomly placed on a circumference
  {
	starting_individual.NewCity(xstart_city,ystart_city,1); //adding the first city          
	double r = 1.; // circumference's radius
	double theta;

	cout << "The cities are randomly placed on a circumference of radius " << r << endl;
	cout << "The position of the first city on the plane is: x = " << xstart_city << " y = " << ystart_city  << endl << endl;

// Generate the starting individual
	for(int icity=2; icity<=ncity; icity++) {
		double x=0.;
		double y=0.;	
		theta = rnd->Rannyu(0., 2.*M_PI);
		x = r*cos(theta);
		y = r*sin(theta);
		starting_individual.NewCity(x, y, icity);
   	}
   }
   else //city randomly placed inside a square
   {
	double l = 2.; // square's side
	cout << "The cities are randomly placed inside a square of side " << l << endl;
// Generate the starting individual
	for(int icity=1; icity<=ncity; icity++) {
		double x=0.;
		double y=0.;	
		x = rnd->Rannyu(-l/2.,l/2.);
		y = rnd->Rannyu(-l/2.,l/2.);
		starting_individual.NewCity(x, y, icity);
   	}
   }

   ReadInput >> peso;

   ReadInput >> prob_crossover;
   ReadInput >> prob_gen_mut1;
   ReadInput >> prob_gen_mut2;
   ReadInput >> prob_gen_mut3;
   ReadInput >> prob_gen_mut4;


   cout << "Exponential selection weight p = " << peso << endl; 
   cout << "Crossover probability = " << prob_crossover*100 << "%" << endl;
   cout << "First genetic mutation probability = " << prob_gen_mut1*100 << "%" << endl;
   cout << "Second genetic mutation probability = " << prob_gen_mut2*100 << "%" << endl;
   cout << "Third genetic mutation probability = " << prob_gen_mut3*100 << "%" << endl;
   cout << "Fourth genetic mutation probability = " << prob_gen_mut4*100 << "%" << endl << endl;

   ReadInput.close();

// generate the first population
   for(int iInd=0; iInd<nindividual; iInd++) {
 	starting_population.NewIndividual(starting_individual);
   }
//shuffle the order of the cities of the starting individuals
   starting_population.StartingPop();

// Check if the starting poulation fulfils the bonds
   if (starting_population.Check()) cout << "Correct initialization" << endl << endl;

// sort the starting population 
   starting_population.SortPop();

//print the starting Population and his best individuals
   if (circle_or_square==1) starting_population.PrintPop("DATA/Circle/Starting_Population.dat");
   if (circle_or_square==0) starting_population.PrintPop("DATA/Square/Starting_Population.dat");
   if (circle_or_square==1) starting_population.PrintBest(0,"DATA/Circle/Best_Individual_Evolution_circle.dat",0);
   if (circle_or_square==0) starting_population.PrintBest(0,"DATA/Square/Best_Individual_Evolution_square.dat",0);
   if (circle_or_square==1) starting_population.PrintHalfAverage("DATA/Circle/Average_Individuals_Evolution_circle.dat",nindividual,0);
   if (circle_or_square==0) starting_population.PrintHalfAverage("DATA/Square/Average_Individuals_Evolution_square.dat",nindividual,0);

//copy the startting population in old population
   for(int iInd=0; iInd<nindividual; iInd++){
   	old_pop.NewIndividual(starting_population.GetIndividual(iInd)); 
   }
}

void Reset(void){

//reset the temp population
   temp.Clear();
//copy old population in temp population
   for(int iInd=0; iInd<nindividual; iInd++){
	temp.NewIndividual(old_pop.GetIndividual(iInd));
   }
}

void Cross_Mut(void){

//individuals selection
   index_imot = temp.Selection(nindividual,rnd,peso);
   index_ifat = temp.Selection(nindividual,rnd,peso);
//crossover
   if ((index_imot!=index_ifat)&&((rnd->Rannyu()) <= prob_crossover))
   {
	temp.CrossOver(rnd->Rannyu(1.,ncity-2.),index_imot,index_ifat);
   }
//first mutation
   if ( rnd->Rannyu() <= prob_gen_mut1)
   {
	temp.Mutation(index_imot,1,ncity,rnd);
	temp.Mutation(index_ifat,1,ncity,rnd);
   }	
//second mutation
  if ( rnd->Rannyu() <= prob_gen_mut2)
   {
	temp.Mutation(index_imot,2,ncity,rnd);
	temp.Mutation(index_ifat,2,ncity,rnd);
   }	
//third mutation
   if ( rnd->Rannyu() <= prob_gen_mut3)
   {
	temp.Mutation(index_imot,3,ncity,rnd);
	temp.Mutation(index_ifat,3,ncity,rnd);
   }	
//forth mutation	
   if ( rnd->Rannyu() <= prob_gen_mut4)
   {
	temp.Mutation(index_imot,4,ncity,rnd);
	temp.Mutation(index_ifat,4,ncity,rnd);
   }
}

void Accumulate(void)
{
// accumulate modified (or not) individuals in new population
   new_pop.NewIndividual(temp.GetIndividual(index_imot));
   new_pop.NewIndividual(temp.GetIndividual(index_ifat));	
}

void PrintResults(int igen)
{
//Print the length of the best path
   if (circle_or_square==1) new_pop.PrintBest(0,"DATA/Circle/Best_Individual_Evolution_circle.dat",igen);
   if (circle_or_square==0) new_pop.PrintBest(0,"DATA/Square/Best_Individual_Evolution_square.dat",igen);

   if(igen%100==0){
	cout << igen << " generation" << endl << endl;
	cout << "----------------------------" << endl << endl;
   }

//Print the averaged length value on the best half of the population
   if (circle_or_square==1) new_pop.PrintHalfAverage("DATA/Circle/Average_Individuals_Evolution_circle.dat",nindividual,igen);
   if (circle_or_square==0) new_pop.PrintHalfAverage("DATA/Square/Average_Individuals_Evolution_square.dat",nindividual,igen);

//Print the best Path at the end of the simulation
   if(igen==ngen){
	if (circle_or_square==1) new_pop.Print_The_Best_Path(0,"DATA/Circle/Best_Individual_circle.dat");
	if (circle_or_square==0) new_pop.Print_The_Best_Path(0,"DATA/Square/Best_Individual_square.dat");
   }
}

void Restart(void)
{
// Reset old population
   old_pop.Clear();
// Copy new population in old population
   for(int iInd=0; iInd<nindividual; iInd++){
	old_pop.NewIndividual(new_pop.GetIndividual(iInd));
   }
// Reset new popolation
   new_pop.Clear();
}




