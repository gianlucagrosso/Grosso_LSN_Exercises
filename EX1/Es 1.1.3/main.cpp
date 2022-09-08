#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


#include "chisquared.h"

using namespace std;
 
int main (int argc, char *argv[]){
   
   if (argc != 3) {
    cout << "Usage: ./main <nchi> <nsubint>  "<< endl;
    return -1;
  }


   int N = 10000; //total numbers to generate

   int M = atoi ( argv[2] ); //number of subintervals

   //file with n° subintervals for data representation

   ofstream file;
   file.open("nsubint.txt");
   file << M << endl; 
   file.close();
   
  //file with 10000 evaluation of chi squared
   
   ofstream file2;
   file2.open("chisquared.txt");
   
   for( int i = 0; i < atoi( argv[1] ); i++ ){
   
      file2 << ChiSquared( M , N) << endl;
   
   }

  
   file2.close();



   return 0;
}