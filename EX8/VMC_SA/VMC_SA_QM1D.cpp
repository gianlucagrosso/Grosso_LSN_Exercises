/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <chrono>
#include <thread> 

#include "VMC_SA_QM1D.h"

using namespace std;
 
int main (int argc, char *argv[]){

  Input(); //VMC Inizialization

  ofstream EnMin("output_EnMin.dat");

  SAReset();
  
  cout << "-------------------------------------------" << endl << endl;

 for(int j = 0; j<10 ; j++){
    temp = T[j];
    beta = 1./temp;
    deltaSA = delta[j];
    SAReset(); //reset SA counters

    cout << "Temperature T: " << T[j] << endl;
    cout << "DeltaSA D: " << deltaSA << endl; 


    for(int i = 0; i < NSA_step; i++){
      
      //Simulated annealing move
      SAMove();
      EnMin << setw(20) << beta << setw(20) << energy << setw(20) << energy_err << setw(20) <<  mu << setw(20) << sigma << setw(20) << SAaccepted/SAattempted << endl; 
     

     
  
      if(j == 9){
        
        if(i == 0) //Save the first value
        {
        en_min = energy;
        en_min_err = energy_err;
        mu_min = mu;
        sigma_min = sigma;
        }

        if(energy <= en_min and ((abs(energy - en_min) >= (en_min_err + energy_err)) or (i == 0))){ //Consider just non compatible energies
        //Save the new min
          en_min = energy;
          en_min_err = energy_err;
          mu_min = mu;
          sigma_min = sigma;

        //Reset input files
          Ham.open("output_Have.dat");
          PsiT_2.open("output_PsiT2.dat");

        //Print config
          Energy(1);
        }
      }

    
    }

    cout << "Acceptance rate " << SAaccepted/SAattempted << endl << endl;
    cout << "-------------------------------------------" << endl << endl;

  }


  cout << "Printing minimal configuration: " << endl;
  cout <<  "Energy: " << en_min << " +/- " << en_min_err << endl;
  cout << "Mu: " << mu_min << " Sigma: " << sigma_min << endl;
  cout << "-------------------------------------------" << endl << endl;



  return 0;
}


void Input(void){


  ifstream ReadInput;


  //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream Seed("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   Seed.close();

  
  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> iQM1D; 
  if(iQM1D == 0){
        cout << "Single quantum particle in 1D      " << endl;
        cout << "Variational Monte Carlo simulation  and Simulated Annealing         " << endl << endl;
        cout << "Trial Wave Function psi_T(x;mu,sigma) = Gauss(mu , sigma) + Gauss(-mu , sigma) " << endl << endl; 
        cout << "External potential V(x) = x^4 - 5/2 * x^2" << endl << endl;
        cout << "The program uses units with hbar = 1 and m = 1  " << endl;
    }

    else{
        cout << "Single quantum particle in 1D      " << endl;
        cout << "Variational Monte Carlo simulation  and Simulated Annealing         " << endl << endl;
        cout << "Trial Wave Function psi_T(x;mu,sigma) = Gauss(mu , sigma) " << endl << endl; 
        cout << "External potential V(x) =  1/2 * x^2" << endl << endl;
        cout << "The program uses units with hbar = 1,  m = 1  and w=1" << endl;
    }

  ReadInput >> mu;
  cout << "The first variational parameter of the trial wave function = " << mu << endl;  
  ReadInput >> sigma;
  cout << "The second variational parameter of the trial wave function = " << sigma << endl;


// starting position
  ReadInput >> start_position;
  x_pos = start_position;

  cout << "The initial position is x = " << x_pos << endl << endl;
  cout << "Initial hamiltonian expectation value = " << DerSecPsi(x_pos , iQM1D) << endl << endl;


  ReadInput >> deltaVMC;

  ReadInput >> nblk;
  Nblk_appo = nblk;

  ReadInput >> nstep;

  ReadInput >> NSA_step;

  cout << "The program perform VMC Metropolis moves with uniform translations" << endl;
  cout << "VMC moves parameter = " << deltaVMC << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  cout << "The program perform SA Metropolis moves with uniform translations" << endl;
  cout << "Number of steps per temperature = " << NSA_step << endl << endl;

//Prepare arrays for measurements
  ih = 0; //Hamiltonian expectation values
  iw = 1; 

  n_props = 2; //Number of observables

  ReadInput >> nbins;				
  ReadInput >> histogram_start;				
  ReadInput >> histogram_end;

// Histogram with the sampled configurations of |psi_trial(x)|^2
  bin_size = (histogram_end - histogram_start)/(double)nbins;

  cout << "The histogram with the sample configurations has " << nbins << " bins " << " with bins dimension = " << bin_size << endl;
  cout << "The histogram starts from x = " << histogram_start << " , and ends to x = " << histogram_end << endl << endl;

  ReadInput.close();

//Reset input files
  Ham.open("output_Have.dat");
  PsiT_2.open("output_PsiT2.dat");

//Evaluate hamiltonian expectation value of the initial position
  VMCMeasure();

//Print initial value for the hamiltonian expectation value
  cout << "Initial hamiltonian expectation value = " << (-1.*DerSecPsi(x_pos , iQM1D)/2. + Pot(x_pos , iQM1D )*Psi(x_pos , iQM1D))/(Psi(x_pos , iQM1D)) << endl << endl;

//Verifing the delta step is the right one
FindDelta();

cout << "The delta step for VMC to have acceptance rate 0.5 is delta = " << deltaVMC << endl << endl;

//Setting up initial temperature and delta 
deltaSA = delta[0];
temp = T[0];
beta = 1./temp;

cout << "Initial SA temperature T = " << temp << endl;
cout << "Initial SA beta B = " << beta << endl;
cout << "Initial deltaSA = " << deltaSA << endl << endl;

//Evaluating initial expectation value
Energy(0);

}




void VMCMove(void)
{ 
  double p;
  double xold,xnew;

//Old
  xold = x_pos;
//New
  xnew = x_pos + deltaVMC*(rnd.Rannyu() - 0.5);

//Metropolis test
  p =  (pow(Psi(xnew, iQM1D),2))/(pow(Psi(xold , iQM1D),2)); 

  if(p>=rnd.Rannyu() or p>=1){
	  x_pos = xnew;
	  VMCaccepted++;
  }
  
  else{
	  x_pos = xold;	
  }

  VMCattempted++;
}


void VMCMeasure()
{
  
  
//reset the hystogram of |psi_trial(x)|^2
  for (int k= 0 ; k<nbins; ++k) walker[iw + k]=0.0;

//update of the histogram of |psi_trial(x)|^2
  for (int i = 0; i<nbins; i++){
	if ( (x_pos > (histogram_start + (i)*bin_size))  && (x_pos < (histogram_start + (i+1)*bin_size))){
		walker[iw +i]++;
	}
  }
  walker[ih] = (-1.*DerSecPsi(x_pos , iQM1D)/2. + Pot(x_pos , iQM1D )*Psi(x_pos , iQM1D))/(Psi(x_pos , iQM1D));
}


void VMCReset(int iblk) //Reset block averages
{ 
  
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props + nbins; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props + nbins; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   VMCattempted = 0;
   VMCaccepted = 0;
}

void VMCAccumulate(void) //Update block averages
{ 

   for(int i=0; i<n_props + nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void VMCAverages(int iblk) //Print results for current block
{   
    double x_current;
    double norma = 0.0; // normalization histogram
    const int wd=12;

    /*cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << VMCaccepted/VMCattempted << endl << endl;*/
    
    Ham.open("output_Have.dat" , ios::app);
    PsiT_2.open("output_PsiT2.dat", ios::app);
    
    stima_ham = blk_av[ih]/blk_norm; // hamiltonian expectation value
    glob_av[ih] += stima_ham;
    glob_av2[ih] += stima_ham*stima_ham;
    err_ham = Error(glob_av[ih],glob_av2[ih],iblk);
    Ham << setw(wd) << iblk <<  setw(wd) << stima_ham << setw(wd) << glob_av[ih]/(double)iblk << setw(wd) << err_ham << endl;


 for (int i=0; i<nbins; i++){
	x_current = histogram_start + (i+0.5)*bin_size;
	stima_PsiT2 = blk_av[iw + i]/blk_norm; // |psi_trial(x)|^2
	glob_av[iw + i] += stima_PsiT2;
	glob_av2[iw + i] += stima_PsiT2*stima_PsiT2;
	err_PsiT2=Error(glob_av[iw + i],glob_av2[iw + i],iblk);
	if (iblk==nblk){
		for (int j=0; j<nbins; j++){
			norma += (glob_av[j + iw]/(double)iblk);
		}
		norma *= bin_size;
		PsiT_2 << setw(wd) << x_current <<  setw(wd) << glob_av[i + iw]/(double)iblk/norma << setw(wd) << err_PsiT2/norma << endl;
	}
    }
    
   //cout << "----------------------------" << endl << endl;

   Ham.close();
   PsiT_2.close();
   
}

void VMCConfFinal(void)
{ 
  ofstream WritePos;

  cout << "Print final position to file pos.final " << endl << endl;
  WritePos.open("pos.final");

  WritePos << " x = " << x_pos << endl;

  WritePos.close();

  rnd.SaveSeed();
}

void Energy(int iPrint){

    x_pos = start_position;

    for(int iblk=1; iblk <= Nblk_appo; ++iblk) //VMC Simulation
	{
		VMCReset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep)
    		{
      			VMCMove();
      			VMCMeasure();
      			VMCAccumulate(); //Update block averages
      		}
    		if(iPrint) VMCAverages(iblk); //Print results for current block
        else SAAverages(iblk);        //Update current energy for SA
        
  	}

}

void FindDelta(void){ //Algorithm to reset deltaVMC in order to always have acceptance rate 0.5
  
  double step = deltaVMC;
 	double rate=0;
	double rate_min=0.47;
	double rate_max=0.52;
	int counter=0;


  do
	{
			counter++;

      VMCaccepted = 0.;
      VMCattempted = 0.;
      x_pos = start_position;

    	for(int istep=1; istep <= 1000; ++istep) VMCMove();
   
        
  	  rate=(double)VMCaccepted/(VMCattempted);

			
			if(rate<rate_min)
			{
				deltaVMC = deltaVMC - deltaVMC/10.; //subtract ten percent
		
			}
			else if(rate>rate_max)
			{
				deltaVMC = deltaVMC + deltaVMC/10.;
			}
			
			if(counter>5000) {
        deltaVMC = step;
        //cout << "Taking too much time to find the right delta step. Exiting!" << endl;
        break;
        }
	}
	while(!(rate_min<rate and rate<rate_max));

}




void SAAverages(int iblk) //Print results for current block
{   
    stima_ham = blk_av[ih]/blk_norm; // hamiltonian expectation value
    glob_av[ih] += stima_ham;
    glob_av2[ih] += stima_ham*stima_ham;
    err_ham = Error(glob_av[ih],glob_av2[ih],iblk);

    if(iblk == Nblk_appo){ 
      energy = glob_av[ih]/(double)iblk;
      energy_err = err_ham;
    }
}

void SAMove(void)
{
  double p; 
  int counter = 0;
  

 //Old parameters
    double energy_old = energy;
    double energy_err_old = energy_err; 
    double mu_old = mu, sigma_old  = sigma;
  
  //New parameters
    mu = mu_old + rnd.Rannyu(-deltaSA , deltaSA);
    sigma = sigma_old + rnd.Rannyu(-deltaSA , deltaSA);

    FindDelta();
    Energy(0);
    Nblk_appo = nblk;

    while(abs(energy - energy_old)< (energy_err + energy_err_old)){ //choose only non compatible energies
      
      Nblk_appo++;
      Energy(0);
      counter++;

      if(counter>200){
        //cout << "Taking too long to find the right energy. Exiting!" << endl;
        break;
      }

    }

//Metropolis test
    p = exp(-beta*(energy - energy_old));


    if(p>=rnd.Rannyu() or p>=1){
	    SAaccepted++;
    }
  
  else{
	  mu  = mu_old;
    sigma = sigma_old;
    energy = energy_old;
    energy_err = energy_err_old; 	
  }

  SAattempted++;


}

void SAReset(void){

  SAaccepted = 0.;
  SAattempted = 0.;

}


double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


double Psi(double x , int iQM1D) //Wave function evaluation 
{   
    
    if(iQM1D) return exp(-pow(x-mu, 2)/(2*pow(sigma, 2)));

    else return exp(-pow(x-mu, 2)/(2*pow(sigma, 2))) + exp(-pow(x+mu, 2)/(2*pow(sigma, 2)));
}

double DerSecPsi(double x, int iQM1D) //Second derivative wave function evaluation 
{   
    if(iQM1D) return exp(-pow(x-mu, 2)/(2*pow(sigma, 2)))*(pow((x-mu)/sigma,2) - 1 )/(pow(sigma, 2));
    
    else return exp(-pow(x-mu, 2)/(2*pow(sigma, 2)))*(pow((x-mu)/sigma,2) - 1 )/(pow(sigma, 2)) + exp(-pow(x+mu, 2)/(2*pow(sigma, 2)))*(pow((x+mu)/sigma,2) - 1 )/(pow(sigma, 2));
}

double Pot(double x , int iQM1D)
{
    if(iQM1D) return pow(x, 2) / 2.;

    else return pow(x ,4) - 5.0*pow(x , 2)/2.0;
}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
