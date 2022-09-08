
#include "Travelling_Salesman_Problem_TSP.h"
#include "funzioni.h"


using namespace std;

int main (int argc, char *argv[]){

// parallelize in the Single Program Multiple Data (SPMD)
   MPI_Init(&argc,&argv); //MPI initialization
   MPI_Comm_size(MPI_COMM_WORLD, &size); //to know how many process(node)
   MPI_Comm_rank(MPI_COMM_WORLD, &node); //to know what node


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
//Exchange best individuals among nodes
	if (igen%nmigration==0) Migration(); //Exchanges among continents	
	Restart(); //Rearrange the populations for the next generational cycle	
   }
			  
   rnd->SaveSeed();
   MPI_Finalize(); //MPI finalization
   return 0;
}


void Input(void){

  ifstream ReadInput;

//Read seed for random numbers
   int p1, p2, temp1, temp2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
	//inizializzo i singoli nodi con diverse coppie di completamento per differenziare le sequenze stocastiche
	for (int irank=0; irank<size; irank++){
		Primes >> temp1 >> temp2;
		if (node==irank){
			p1=temp1;
			p2=temp2;
		}
	}
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

   if (node==0)
   {
	cout << endl;
	cout << "This program solve the Travelling salesman Problem with a parallelized Genetical Algorithm" << endl;
	cout << "Each node perform an independent GA search: the so-called 'Continents' " << endl;
	cout << "The bonds are: the salesman must visit one and only one every city and must be back to the fisrt city at the end of the path!" << endl << endl;
   }


  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> ncity;
  if (node==0) cout << "The path is made of " << ncity << " cities" << endl;

  ReadInput >> nindividual;
  if (node==0) cout << "Each population is made of " << nindividual << " individuals" << endl <<endl;

  ReadInput >> ngen;
  if (node==0) cout << "The generations (new populations) produced by tha algorithm are : " << ngen << endl;

  ReadInput >> nmigration;
  if (node==0) cout << "Every " << nmigration << " generations the Continents exchange their best individuals randomly" << endl;

  ReadInput >> xstart_city; 
  ReadInput >> ystart_city;

  ReadInput >> circle_or_square;

//Generate the starting individual for the processors marked with node=0
  if (node==0)
  {    
      if(circle_or_square==0) //city randomly placed inside a square
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

      if(circle_or_square == 2) //American capital cities
      {  
         ifstream Capitals("American_capitals.dat");

         double x , y;
         int icity = 1; 

         while(!(Capitals.eof())){
            Capitals >> x >> y; 
            starting_individual.NewCity(x, y, icity++);
            }
      }



   }

   ReadInput >> peso;

   ReadInput >> prob_crossover;
   ReadInput >> prob_gen_mut1;
   ReadInput >> prob_gen_mut2;
   ReadInput >> prob_gen_mut3;
   ReadInput >> prob_gen_mut4;

   if (node==0) cout << "Il peso esponenziale nella selezione p = " << peso << endl;
   if (node==0) cout << "Crossover probability = " << prob_crossover*100 << "%" << endl;
   if (node==0) cout << "First genetic mutation probability = " << prob_gen_mut1*100 << "%" << endl;
   if (node==0) cout << "Second genetic mutation probability = " << prob_gen_mut2*100 << "%" << endl;
   if (node==0) cout << "Third genetic mutation probability = " << prob_gen_mut3*100 << "%" << endl;
   if (node==0) cout << "Fourth genetic mutation probability = " << prob_gen_mut4*100 << "%" << endl << endl;

   ReadInput.close();

//Broadcasts the generated starting individual from the root 0 to all the other processes in the default communicator
   double x_starting_city[ncity];
   double y_starting_city[ncity]; 
   int index_starting_city[ncity]; 

   for(int icity=0; icity<ncity; icity++){
	if (node == 0){
		x_starting_city[icity] = starting_individual.GetCity(icity).GetX();
		y_starting_city[icity] = starting_individual.GetCity(icity).GetY();
		index_starting_city[icity] = starting_individual.GetCity(icity).GetIndex();
	}
	else{
		x_starting_city[icity] = 0.;
		y_starting_city[icity] = 0.;
		index_starting_city[icity] = 0;
	}		
   }

   MPI_Bcast(x_starting_city,ncity,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(y_starting_city,ncity,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(index_starting_city,ncity,MPI_INT,0,MPI_COMM_WORLD);

   for (int irank=1; irank<size; irank++){
	if (node==irank){
		for(int icity=0; icity<ncity; icity++) {
			starting_individual.NewCity(x_starting_city[icity],y_starting_city[icity],index_starting_city[icity]);
   		}
	}
   }

//Generate the first population
   for(int iInd=0; iInd<nindividual; iInd++) {
 	starting_population.NewIndividual(starting_individual);
   }
//Shuffle the order of the cities of the starting individuals in each processors
   for (int irank=0; irank<size; irank++){
	int shuffle_rip = 0;
	if (node==irank){
		while(shuffle_rip<irank+1){ // shuffle_rip = number of shuffle  
			starting_population.StartingPop();
			shuffle_rip++;
		}
	}
   }

// Check if the starting poulation fulfils the bonds
   if (starting_population.Check()){
	if (node==0) cout << "Correct initialization of every node" << endl << endl;
   }

// sort the starting population 
   starting_population.SortPop();

//Printing path
if(circle_or_square==2)   path = "DATA/AmericanCities/Node " + to_string(node) + "/";
if(circle_or_square==1)   path = "DATA/Circle/Node " + to_string(node) + "/";
if(circle_or_square==0)   path = "DATA/Square/Node " + to_string(node) + "/";

//   if (node==0) cout << "The first random population: " << endl;
//..starting population
starting_population.PrintPop((path + "Starting_Population.dat").c_str());
//...the length of the best path of starting population for each processor 
starting_population.PrintBest(0,(path+"Best_Individual_Evolution.dat").c_str(),0);
//...the averaged length value on the best half of the population for each processor 
starting_population.PrintHalfAverage((path + "Average_Individuals_Evolution.dat").c_str(),nindividual,0);
   

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
//Print the length of the best path for each processor
   new_pop.PrintBest(0,(path + "Best_Individual_Evolution.dat").c_str(),igen);

   if (node==0){
   	if(igen%100==0){
		cout << igen << " generation" << endl;
		cout << "----------------------------" << endl << endl;
   	}
   }
//Print the averaged length value on the best half of the population
  new_pop.PrintHalfAverage((path + "Average_Individuals_Evolution.dat").c_str(),nindividual,igen);

//Print the best Path at the end of the simulation
   if(igen==ngen){
	   new_pop.Print_The_Best_Path(0,(path + "Best_Individual.dat").c_str());
   }

//Printing best of the best 
   vector<double> length(size);
   double ilength = new_pop.GetIndividual(0).CostFun();

   MPI_Gather(&ilength , 1 , MPI_DOUBLE , &length[node] , 1 , MPI_DOUBLE , 0, MPI_COMM_WORLD );

   if(node==0){

      ofstream Best_of_the_Best;

      if(circle_or_square==0) Best_of_the_Best.open("DATA/Square/Best_of_the_Best.dat", ios::app);
      if(circle_or_square==1) Best_of_the_Best.open("DATA/Circle/Best_of_the_Best.dat", ios::app);
      if(circle_or_square==2) Best_of_the_Best.open("DATA/AmericanCities/Best_of_the_Best.dat", ios::app);  

      int bestnode = min_element(length.begin() , length.end()) - length.begin();

      Best_of_the_Best << setw(12) << igen << setw(12) <<  *(min_element(length.begin() , length.end())) << setw(12) << bestnode << endl;

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

void Migration(void){

   MPI_Status stat1, stat2, stat3, stat4, stat5, stat6;

//contenitori ausiliari
   double x_city_1[ncity]; 
   double y_city_1[ncity]; 
   int index_city_1[ncity]; 

   double x_city_2[ncity]; 
   double y_city_2[ncity]; 
   int index_city_2[ncity];  

//max exchanges among continents: size
 for (int j=0; j<size; j++)
 {
//randomly selection of the rank
   for(int i=0; i<2; i++){
	if (node==0){   
		rk[i] = (int)(size*rnd->Rannyu());
		//cout << "Rank " << rk[i] << endl;
   	}
	else{
		rk[i] = 0;
	}
   }

   MPI_Bcast(rk,2,MPI_INT,0,MPI_COMM_WORLD);
   
 
// Send and Receive
   if (rk[0] != rk[1] )
   {
	for(int icity=0; icity<ncity; icity++){
		if (node==rk[0]){
			x_city_1[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetX();
			y_city_1[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetY();
			index_city_1[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetIndex();
		}
   		if (node==rk[1]){
			x_city_2[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetX();
			y_city_2[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetY();
			index_city_2[icity] = ((new_pop.GetIndividual(0)).GetCity(icity)).GetIndex();
		}
	}
// x coordinate cities
	if (node==rk[0]){
		MPI_Send(&x_city_1[0],ncity,MPI_DOUBLE,rk[1],itag1,MPI_COMM_WORLD);
		MPI_Recv(&x_city_2[0],ncity,MPI_DOUBLE,rk[1],itag2,MPI_COMM_WORLD,&stat2);
   	}
   	if (node==rk[1]){
		MPI_Recv(&x_city_1[0],ncity,MPI_DOUBLE,rk[0],itag1,MPI_COMM_WORLD,&stat1);
		MPI_Send(&x_city_2[0],ncity,MPI_DOUBLE,rk[0],itag2,MPI_COMM_WORLD);
   	}
// y coordinate cities
	if (node==rk[0]){
		MPI_Send(&y_city_1[0],ncity,MPI_DOUBLE,rk[1],itag3,MPI_COMM_WORLD);
		MPI_Recv(&y_city_2[0],ncity,MPI_DOUBLE,rk[1],itag4,MPI_COMM_WORLD,&stat4);
	}
	if (node==rk[1]){
		MPI_Recv(&y_city_1[0],ncity,MPI_DOUBLE,rk[0],itag3,MPI_COMM_WORLD,&stat3);
		MPI_Send(&y_city_2[0],ncity,MPI_DOUBLE,rk[0],itag4,MPI_COMM_WORLD);
	}
// index cities
	if (node==rk[0]){
		MPI_Send(&index_city_1[0],ncity,MPI_INT,rk[1],itag5,MPI_COMM_WORLD);
		MPI_Recv(&index_city_2[0],ncity,MPI_INT,rk[1],itag6,MPI_COMM_WORLD,&stat6);
	}
	if (node==rk[1]){
		MPI_Recv(&index_city_1[0],ncity,MPI_INT,rk[0],itag5,MPI_COMM_WORLD,&stat5);
		MPI_Send(&index_city_2[0],ncity,MPI_INT,rk[0],itag6,MPI_COMM_WORLD);
	} 

	for(int icity=0; icity<ncity; icity++){
		if (node==rk[0]){
			new_pop.Setcity(0,icity,x_city_2[icity],y_city_2[icity],index_city_2[icity]);
		}
		if (node==rk[1]){
			new_pop.Setcity(0,icity,x_city_1[icity],y_city_1[icity],index_city_1[icity]);
		}
	}
   }
 }   	
}


