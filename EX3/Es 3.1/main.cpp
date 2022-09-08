#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <functional>

#include "random.h"
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){

if (argc != 3) {
        cout << "Usage: ./main <nsteps> <nblocks>"<< endl;
        return -1;
    }


//setting up variables
    int Nsteps = atoi ( argv[1] ); //number of total values 
    int Nblocks = atoi ( argv[2] ); //number of blocks for Data Blocking




    double t=0;  //time of the sell of the option
    double S0=100;  //initial asset prize
    double T=1;    //delivery time
    double K=100;   //strike price
    double r=0.1;  //drift and risk-free constant rate
    double sigma=0.25;  //volatility

    Random rnd;

    vector <double> C_appo(Nsteps , 0.);
    vector <double> P_appo(Nsteps , 0.);

    vector <double> Call(Nblocks); //running averages for call prices
    vector <double> Call_error(Nblocks);//running errors for call prices

    vector <double> Put(Nblocks);   //analogous for put options
    vector <double> Put_error(Nblocks);

    double S = 0.;  //price at time T


////DIRECT SAMPLING


    for (int i = 0; i < Nsteps; i++){

			S=S0*exp((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,1)*sqrt(T)); //sample asset price at final T 
			
            if((S-K)>0){//check wheter the holder uses the option or not, both in call and put case
				
                C_appo[i] = exp(-r*T)*( S - K);
                P_appo[i] = 0.;
            }

			else{

				P_appo[i] = exp(-r*T)*(K - S);
                C_appo[i] = 0.; 
            
            }

		}

    Call = Average(C_appo, Nblocks);
    Call_error = DevStd(C_appo , Nblocks);

    Put = Average(P_appo , Nblocks);
    Put_error = DevStd(P_appo , Nblocks);

    Print(Call.begin() , Call.end() , "DATA/call_direct.txt");
    Print(Call_error.begin() , Call_error.end() , "DATA/dev_call_direct.txt");

    Print(Put.begin() , Put.end(), "DATA/put_direct.txt" );
    Print(Put_error.begin() , Put_error.end(), "DATA/dev_put_direct.txt" );
    



////DISCRETE TIME EVOLUTION SAMPLING

    double passo=(T-t)/100.;

    for (int i = 0; i < Nsteps; i++){

            t = 0.;
            S = S0;

			while (t<T) //evolve the asset price step by step
			{
				t+=passo;
				S=S*exp((r-0.5*sigma*sigma)*passo+sigma*rnd.Gauss(0,1)*sqrt(passo));
			}
			
            if((S-K)>0){//check wheter the holder uses the option or not, both in call and put case
				
                C_appo[i] = exp(-r*T)*( S - K);
                P_appo[i] = 0.;
            }

			else{

				P_appo[i] = exp(-r*T)*(K - S);
                C_appo[i] = 0.; 
            
            }

	}


    Call = Average(C_appo, Nblocks);
    Call_error = DevStd(C_appo , Nblocks);

    Put = Average(P_appo , Nblocks);
    Put_error = DevStd(P_appo , Nblocks);

    Print(Call.begin() , Call.end() , "DATA/call_discrete.txt");
    Print(Call_error.begin() , Call_error.end() , "DATA/dev_call_discrete.txt");

    Print(Put.begin() , Put.end(), "DATA/put_discrete.txt" );
    Print(Put_error.begin() , Put_error.end(), "DATA/dev_put_discrete.txt" );
    

    rnd.SaveSeed();

    return 0;


}