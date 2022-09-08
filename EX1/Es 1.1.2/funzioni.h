#pragma once

#include <iostream> 
#include <vector>    
#include <algorithm>   
#include <iterator>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "random.h"

using namespace std; 

//uniform number's vector generator 

template <typename T> vector<T> VectUnif( int M ){

  Random rnd;


  vector<T> r(M, 0.);

  for(int i=0; i<M; i++){
    
    if (typeid(T) == typeid(double)) r[i] = rnd.Rannyu();
    else r[i] = (int)(rnd.Rannyu() * INT64_MAX);
   
  }


  rnd.SaveSeed();


  return r; 


}

template <typename T> vector<T> VectUnif( int M , T a, T b ){

  Random rnd;


  vector<T> r(M, 0.);

  for(int i=0; i<M; i++){
    
    if (typeid(T) == typeid(double)) r[i] = rnd.Rannyu(a, b);
    else r[i] = (int)(rnd.Rannyu(a, b));
   
  }


  rnd.SaveSeed();


  return r; 


}


//data blocking average and dev.std

template <typename T> vector<T> Average(  const vector<T>& r , int N ) {

      
  vector<T> ave(N, 0.);
    
  int L=static_cast<int>(r.size())/N;  //lenght of block for data blocking

  for(int i=0; i<N; i++){
    
      T sum1 = 0; 
             
      for(int j=0; j<(i+1)*L; j++){
           
            sum1 += r[j];   //sum r_j {j=0, end of block i}
            
            }
        
      ave[i] = sum1/((i+1)*L); //cumulative average
                
              
      }

    return ave; 
       
};



template <typename T> vector<T> DevStd(  const vector<T>& r, int N ) {

      
  vector<T> ave(N, 0.);
  vector<T> ave2(N, 0.);
  vector<T> sum_prog(N, 0.);
  vector<T> sum2_prog(N, 0.); 
  vector<T> err_prog(N, 0.);

  int L=static_cast<int>(r.size())/N;  //lenght of block for data blocking
  

  for(int i=0; i<N; i++){
      T sum1 = 0; 
             
      for(int j=0; j<L; j++){

          int k = j+i*L;   //sum of variables in each block
          sum1 += r[k];

            }
        
      ave[i] = sum1/L;    //mean in each block
      ave2[i] = (ave[i])*(ave[i]); //mean squared in each block
        

      for(int j = 0; j < i+1; j++ ){

            sum_prog[i] += ave[j];  //sum of A_j {0, block n° i}
            sum2_prog[i] += ave2[j]; // sum of A_j^2  ""

               }

      sum_prog[i] /= (i+1); //cumulative average
      sum2_prog[i] /= (i+1);  //cumulative square average

      if( i != 0) err_prog[i] = sqrt( (sum2_prog[i]-sum_prog[i]*sum_prog[i])/i );  //progressive SDM (with N-1 instead of N, bc we have a relative small number of blocks)
      else err_prog[i] = 0.; //first dev. std is 0 since we have just a single value 
        
       }
    
   return err_prog;

};



//screen printing of vector's content

template <typename T> void Print( T __start , T __end ) {

	int ncol;	
  cout << endl << "Type the number of columns desired for printing on screen: " << endl;
  cin >> ncol;

  int appo = int(__end - __start) % ncol; //number of complete lines 

  cout << "Dataset: " << endl;
  	
  for(auto it =__start ; it!=__end - appo; it = it + ncol){
	  
    for(auto it2 = it ; it2!=it + ncol; it2++)
  		cout << setw(12) << *it2;
	
  cout << endl;
	
  }
	
  for(auto it = __end - appo; it!=__end; it++ )  //printing the residual stuff
  	cout << setw(12) << *it;
  
 

  cout<<endl;
 

};


//overloading of previous function to print on file

template <typename T> void Print( T _start , T _end , const char * filename) {
  
  ofstream outputfile (filename);
    
  if ( !outputfile ) {
    cerr << "Can't create file" << filename << endl;
    exit(2);
  }

int ncol = 1;
 /* cout << endl <<  "Type the number of columns desired for printing on file: " << endl;
  cin >> ncol; */

int appo = int(_end- _start) % ncol; 

  
for(auto it =_start ; it!=_end - appo; it = it + ncol){
	for(auto it2 = it ; it2!=it + ncol; it2++)
  		outputfile << setw(8) << *it2;
	outputfile << endl;
	}
	
  for(auto it = _end - appo; it != _end; it++ )
  	outputfile << setw(8) << *it;

  outputfile.close();
  
  };


