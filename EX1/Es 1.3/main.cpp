
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


#include "random.h"
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){
   
    if (argc != 3) {
    cout << "Usage: ./main <nblocks> <nthrows>  "<< endl;
    return -1;
  }

    Random rnd; 


    //saving the number of blocks and throws per block for data representation 
    int nblocks = atoi(argv[1]);
    int nthrows = atoi(argv[2]);

    ofstream file;
    file.open("ndata.txt");
    file << nblocks<< endl << nthrows << endl; 
    file.close();

    //lenght of the needle (L) and space btw the lines (d)
    double L = 3.;
    double d = 5.; 

    //centre of the needle, orientation, y of the edges
    double y_centre , sinT , y1 , y2 ;

    vector<double> pi ( nblocks , 0. );
    vector<double> devstd (nblocks, 0. );

    int Nhit = 0; 


    //measure: the hit is counted if the edges are in the two different half planes defined by the line y=0 in the middle
    for (int i = 0; i < nblocks; i++ ){

        Nhit = 0; 

        for(int j = 0; j < nthrows; j++ ){

            y_centre = rnd.Rannyu(-d/2. , d/2.);
            sinT = sin(rnd.Unif_Theta());
            y1 = y_centre + (L/2.)*sinT;
            y2 = y_centre - (L/2.)*sinT;
                           
            if((y1>=0 && y2<=0) || (y1<=0 && y2>=0)) Nhit++; 
          
            
        }

        pi[i] = (double)(2. * L * nthrows)/(double)(Nhit* d);

        
    }

    

    //data blocking
    devstd = DevStd(pi , nblocks);
    pi = Average( pi , nblocks);

    Print( pi.begin(), pi.end(), "pi.txt"); 
    Print ( devstd.begin(), devstd.end(), "error.txt"); 

    rnd.SaveSeed();
    
return 0;

}