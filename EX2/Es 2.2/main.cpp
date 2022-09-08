#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <functional>

#include "random.h"
#include "funzioni.h"
#include "RWalk.h"

using namespace std;



int main (int argc, char *argv[]){


//verifying the discrete probability distribution for discrete RW
//each event has a third as probability of occurence

    vector<double> p(6 , static_cast<double>(1)/6);
    vector<double> v (10000 , 0.);
    
    Random rnd; 

    for (int i = 0; i < static_cast<int>(v.size()); i++)
        v[i] = rnd.Discrete(p);

    
    Print(v.begin() , v.end() , "DATA/discrete.txt");
    
    rnd.SaveSeed();
    

//setting up variables 
    if (argc != 4) {
        cout << "Usage: ./main <nRW> <nsteps> <nblocks>"<< endl;
        return -1;
    }

    int NRWalk = atoi( argv[1] ); //number of Random Walks

    int Nsteps = atoi ( argv[2] ); //number of steps for each RW

    int Nblocks = atoi ( argv[2] ); //number of blocks for Data Blocking

    vector<double> O(3, 0.);  //orgin point of the Rwalk

    RWalk rw( 1. , 3 , O);

////DISCRETE 3D RANDOM WALK

    
//support variables
    vector<vector<double>> appo(NRWalk);
    vector<double> appo2 (NRWalk);

   
    vector<double> discrete_distance (Nsteps);
    vector<double> dev_discrete_distance (Nsteps);
    
//moving each NRWalk-th Random Walk of one step and saving the distance
//appo is a vector of the distances from the origin of NRwalk Random Walk of Nsteps steps

    for ( int j = 0; j < NRWalk; j++){
        
        for ( int i =0; i < Nsteps; i++ ){
        
            rw.RW_move_discrete();
            appo[j].push_back( rw.get_norm() );

        }

        rw.move_to_origin();
    }

//saving the distance ad each step for every RW in appo2 and performing data blocking
    for(int l =0; l < Nsteps; l++){

        for( int k = 0; k < NRWalk; k++ )  

            appo2[k] = appo[k][l];


        discrete_distance[l] = sqrt(Average(appo2 , Nblocks)[Nblocks - 1]); //taking the last cumulative average and devstd in data blocking
        dev_discrete_distance[l] = DevStd (appo2 , Nblocks)[Nblocks - 1] * 0.5 / discrete_distance[l]; //error propagation err(sqrt(x))=0.5/sqrt(x)
    }

    Print(discrete_distance.begin() , discrete_distance.end() , "DATA/discreteRW.txt" );
    Print(dev_discrete_distance.begin() , dev_discrete_distance.end() , "DATA/dev_discreteRW.txt" );
    



////CONTINUOS RANDOM WALK


//support variables
    vector<vector<double>> appo3(NRWalk);
    vector<double> appo4 (NRWalk);

   
    vector<double> continuos_distance (Nsteps);
    vector<double> dev_continuos_distance (Nsteps);
    
//moving each NRWalk-th Random Walk of one step and saving the distance
//appo is a vector of the distances from the origin of NRwalk Random Walk of Nsteps steps

    rw.move_to_origin();

    for ( int j = 0; j < NRWalk; j++){
        
        for ( int i =0; i < Nsteps; i++ ){
        
            rw.RW_move_continuos();
            appo3[j].push_back( rw.get_norm() ) ;
        }

        rw.move_to_origin();
    }

//saving the distance ad each step for every RW in appo2 and performing data blocking
    for(int l =0; l < Nsteps; l++){

        for( int k = 0; k < NRWalk; k++ )   

            appo4[k] = appo3[k][l];

        continuos_distance[l] = sqrt(Average(appo4 , Nblocks)[Nblocks - 1]); //taking the last cumulative average and devstd in data blocking
        dev_continuos_distance[l] = DevStd(appo4 , Nblocks)[Nblocks - 1] * 0.5 / continuos_distance[l]; //error propagation
    }

    Print(continuos_distance.begin() , continuos_distance.end() , "DATA/continuosRW.txt" );
    Print(dev_continuos_distance.begin() , dev_continuos_distance.end() , "DATA/dev_continuosRW.txt" );
    
    
    
    return 0;





}