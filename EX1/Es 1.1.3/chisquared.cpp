
#include "chisquared.h"



double ChiSquared( int nsubint , int nthrows){

    vector<int> counts (nsubint , 0); 
    vector<double> appo (nsubint , 0); 
    double chisq = 0.;

    Random rnd; 


    for(int i = 0; i < nthrows; i++ )
        counts[rnd.Rannyu()*nsubint]+=1;

    //Print(counts.begin() , counts.end());

    

    transform(counts.begin(), counts.end(), appo.begin(), [&]( int Oi ){  return (double)pow(Oi-nthrows/nsubint , 2)/(nthrows/nsubint); });
    
    
    chisq=accumulate(appo.begin(), appo.end() , 0.0);

    rnd.SaveSeed();


    return chisq;

};

