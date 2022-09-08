
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){
   
   if (argc != 2) {
    cout << "Usage: ./main <nSum>  "<< endl;
    return -1;
  }

    int nsum = atoi(argv[1]); //number of the sums

    int N[] = {1 , 2 , 10 , 100 }; //number of variables in a sum

//file with n° sums for data representation

    ofstream file;
    file.open("nsum.txt");
    file << nsum << endl; 
    file.close();

//files with the data for the histos 

    Random rnd; 

    ofstream file1, file2, file3;
    file1.open("UniformDistr.txt");
    file2.open("ExpDistr.txt");
    file3.open("CauchyDistr.txt");

    double  SUn = 0. ,  SExp = 0. , SCau = 0. ; //appo variables

    for (int i = 0; i < 4; i++){

        
    
        for(int j = 0; j < nsum; j++ ){

            SUn = SExp = SCau = 0.;
        
            for( int l = 0; l < N[i]; l++ ){ //looping on N to get sums with different number of addends to test the central limit teo. 
 
                 
                SUn += rnd.Rannyu();
                SExp += rnd.Exp(1.); 
                SCau += rnd.Cauchy(0. , 1.);  

            }

            file1 << SUn/N[i] << endl;
            file2 << SExp/N[i] << endl;
            file3 << SCau/N[i] << endl; 
        
        }

    }

    file1.close();
    file2.close();
    file3.close();

return 0; 

}