#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <functional>

#include "random.h"
#include "funzioni.h"

using namespace std;

double pi = M_PI;

double integrand(double x) { return (pi/2)*cos((pi*x)/2);} 

//function to implement the antithetic variates method --> its integral in [0,1] is the same as the function above

double integrand2(double x) { return (pi/4)*( cos((pi*x)/2) + cos((pi*(1-x))/2) );} 

//p.d.f  for importance sampling

double myfun(double x) { return (1- pi * pi / 8.* pow(x - 0.5 , 2)) / (1 - pi * pi / 96.);  }

//double myfun2(double x){ return 3*(1-x*x)/2; }

 





int main (int argc, char *argv[]){

    if (argc != 3) {
        cout << "Usage: ./main <nsteps> <nblocks>"<< endl;
        return -1;
    }



//setting up the random generator

    Random rnd;

//defining variables

    int Nsteps = atoi( argv[1] ); //number of steps

    int Nblocks = atoi ( argv[2] ); //number of blocks

//saving the number of blocks and throws per block for data representation 

    ofstream file;
    file.open("DATA/ndata.txt");
    file << Nblocks<< endl << Nsteps << endl; 
    file.close();

//defining variables

    vector<double> I(Nblocks);
    vector<double> sigmaI(Nblocks);

    vector<double> appo(Nsteps , 0. );





////UNIFORM SAMPLING
//generating Nsteps value of integrand for the uniform method

    for ( int i = 0; i < Nsteps; i++ )
        appo[i] = integrand(rnd.Rannyu()); 


//evaluating the integral through data blocking

    I = Average<double>(appo , Nblocks); 
    sigmaI = DevStd<double>(appo , Nblocks);

    Print(I.begin() , I.end() , "DATA/unif_integral.txt");
    Print(sigmaI.begin() , sigmaI.end() , "DATA/dev_unif_integral.txt");




////IMPORTANCE SAMPLING + ANTITHETIC VARIATES METHOD
//generating Nsteps value of integrand2/myfun where myfun is the pdf
//used for importance sampling
    
    double xi = 0.;
        
    for ( int i = 0; i < Nsteps; i++ ){
        
        //sampling of the pdf through Accept / Reject method
        xi = rnd.AR(0. , 1. ,  myfun(0.5) , myfun ); //max value of myfun is in x=1/2
        appo[i] = integrand2(xi)/myfun(xi); 
    

    }

    I = Average<double>(appo , Nblocks); 
    sigmaI = DevStd<double>(appo , Nblocks);

    Print(I.begin() , I.end() , "DATA/IS_integral.txt");
    Print(sigmaI.begin() , sigmaI.end() , "DATA/dev_IS_integral.txt");
    

    //saving the random variable to verify its pdf is the one expected
    
    vector<double> x ( pow(10 , 5) , 0.);

    for (int i = 0; i < pow(10 , 5); i++ )
        x[i] = rnd.AR(0. , 1. ,  myfun(0.5) , myfun );
    
    Print(x.begin() , x.end() , "DATA/AR_pdf.txt");

    rnd.SaveSeed();

   return 0;
   
   }

