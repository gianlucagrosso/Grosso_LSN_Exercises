#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){
   
   if (argc != 3) {
    cout << "Usage: ./main <nsteps> <nblocks>"<< endl;
    return -1;
  }


   int M = atoi( argv[1] ); //number of steps

   int N = atoi ( argv[2] ); //number of blocks

   
   //file with nsteps and nblocks for data representation
   
   ofstream file;
   file.open("ndata.txt");
   file << M << endl << N << endl; 
   file.close();


   
   vector<double> r = VectUnif<double>(M); 
   
   vector<double> ave = Average<double> (r , N);
   vector<double> err = DevStd<double> (r, N); 

   
    
   Print( ave.begin(), ave.end(), "average.txt"); 
   Print ( err.begin(), err.end(), "error.txt"); 

  

  
   return 0;
}