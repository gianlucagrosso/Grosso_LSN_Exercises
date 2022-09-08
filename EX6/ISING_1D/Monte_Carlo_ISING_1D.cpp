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
#include <ostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization

  if(measure==1 and restart==0){ //Equilibration 
    for(int i = 0; i < 50 * nequi; i++)
      Move(metro);
    }


  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  

  return 0;
}


void Input(void)
{
  ifstream ReadInput , Seed , FinalConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;

  std::ostringstream appostr;
  appostr << std::fixed << setprecision(1) << temp;
  tempstr = appostr.str();

  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> restart;

  ReadInput >> measure;

  ReadInput >> nequi;

  ReadInput >> nblk;

  ReadInput >> nstep;
  
  if(metro==1){ 
    cout << "The program perform Metropolis moves" << endl;
    alg = "metro";
    path = "DATA/"+alg+"/";
  }

  else{ 
    cout << "The program perform Gibbs moves" << endl;
    alg = "gibbs";
    path = "DATA/"+alg+"/";
  }

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   if(restart) Seed.open("seed.out");
   else Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   Seed.close();




//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration

if(restart)
{
  FinalConf.open(path + "config_T= " + tempstr + ".final");
  for (int i=0; i<nspin; i++)
    FinalConf>>s[i];
}

else{
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial energy           = " << walker[iu]/(double)nspin << endl;
  cout << "Initial magnetization  = " << walker[im]/(double)nspin << endl;
  cout << "Initial magnetic susceptibility     = " << walker[ix]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted++;

    if(metro==1) //Metropolis
    {
      sm = -1.*s[o]; //flip the selected spin
      energy_old = Boltzmann(s[o] , o); //old config energy
      energy_new = Boltzmann(sm , o);  //new config energy
      p = exp(-beta* (energy_new - energy_old)); //ratio of boltzann weights
      
      if(p>1 or rnd.Rannyu() <= p){ //accept only if the probability is higher than 1 or if we extract a number below the prob. threshold
        
        s[o] *= -1;
        accepted++;

      }
      else ;
    
    }
    
    else //Gibbs sampling
    {
      energy_up = Boltzmann(+1 , o); //s_k=+1 config energy
      energy_down = Boltzmann(-1 , o);  //s_k=-1 config energy
      p = 1/(1+exp(-beta*(energy_up - energy_down))); //Prob(s[o]=-1|...)=p and Prob(s[o]=+1|...)=1-p
      
      if(rnd.Rannyu()>=p){
        
        s[o] = 1;
        accepted++;

      }

      else{

        s[o] = -1;
        accepted++;

      }

    }


  }
}

double Boltzmann(int sm, int ip)
{ //the energy change is caused by the flipped spin sm only
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i]; 
  }

  walker[iu] = u;
  walker[ic] = u*u;

  walker[im] = m; 
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

  if(h == 0.0)
  {  
    if(measure == 1)  Ene.open(path + "measure/output.ene.dat",ios::app);
    else  Ene.open(path + "equilibration/output.ene_T=" + tempstr + ".dat",ios::app);
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    
    if(measure == 1 and iblk == nblk) Ene <<  setw(wd) << temp << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    else if(measure == 0) Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    
    Ene.close();



    if(measure == 1)  Heat.open(path + "measure/output.heat.dat",ios::app);
    else  Heat.open(path + "equilibration/output.heat_T=" + tempstr + ".dat",ios::app);

    stima_c = pow(beta, 2) * (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; //heat capacity per spin
    glob_av[ic]  += stima_c;
	  glob_av2[ic] += stima_c*stima_c;
	  err_c = Error(glob_av[ic],glob_av2[ic],iblk);

    if(measure == 1 and iblk == nblk) Heat <<  setw(wd) << temp << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    else if(measure == 0) Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  	
    Heat.close();



    if(measure == 1)  Chi.open(path + "measure/output.chi.dat",ios::app);
    else  Chi.open(path + "equilibration/output.chi_T=" + tempstr + ".dat",ios::app);
		
    stima_x = (beta)*blk_av[ix]/blk_norm/(double)nspin; //magnetic susceptibility
		glob_av[ix] += stima_x;
		glob_av2[ix] += stima_x*stima_x;
		err_x=Error(glob_av[ix], glob_av2[ix],iblk);

    if(measure == 1 and iblk == nblk) Chi <<  setw(wd) << temp << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    else if(measure == 0) Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  	
    Chi.close();

}

else
{
  if(measure == 1) Mag.open(path + "measure/output.mag.dat",ios::app);  //magnetization
  else Mag.open(path + "equilibration/output.mag_T=" + tempstr + ".dat",ios::app);	

  stima_m = blk_av[im]/blk_norm/(double)nspin; 
  glob_av[im] += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
  
  if(measure == 1 and iblk == nblk) Mag << setw(wd) << temp << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  else if(measure == 0) Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  
  Mag.close();

}

    cout << "----------------------------" << endl << endl;

}



void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(path + "config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
