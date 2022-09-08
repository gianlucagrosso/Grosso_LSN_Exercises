/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"

using namespace std;

Random :: Random(){

   int seed[4];
   int p1, p2;

   ifstream Primes("Primes");

   if (Primes.is_open()){
   
      Primes >> p1 >> p2 ;
   
   } else cerr << "PROBLEM: Unable to open Primes" << endl;

   Primes.close();

   ifstream input("seed.in");
   ifstream input2 ("seed.out"); 

   string property;


   if(input2){

      while ( !input2.eof() ){
        
            input2 >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
         
      }
    
      input2.close();

      }

   else{
      if (input.is_open()){

      while ( !input.eof() ){
        
        input >> property;
        
        if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
         }
      }
    
      input.close();
   
      }else cerr << "PROBLEM: Unable to open seed.in" << endl;

   }

   
}

Random :: ~Random(){}

void Random :: ResetSeed(){
   remove("seed.out");
}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double decayrate) {
   return -log(1-Rannyu())/decayrate; 
}

double Random :: Cauchy(double mean, double width){
   return width*tan(M_PI*(Rannyu()-0.5)) + mean; 
}

double Random:: AR(double a , double b , double pmax , function<double (double)> p ){
   
   double x = 0.;
   
   do{
      
      x = Rannyu(a , b);

   }while(Rannyu() >= (p(x)/pmax));

   //cout << x << endl;

   return x; 

}

int Random :: Discrete ( const vector<double>& p ){

   double x = Rannyu();

   
   function<int(int , int)> binarySearch = [&]( int low , int high ) -> int{
      
      if ( high >= low ){

         int mid = low + (high - low) / 2;
         
         if ( x > accumulate(p.begin() , next(p.begin() , ( mid - 1 ) + 1 ) , 0.0 ) and x <= accumulate(p.begin() , next(p.begin() , mid + 1) , 0.0 ) )
               return mid; 

         if ( x > accumulate(p.begin() , next(p.begin() , mid  + 1) , 0.0 ) )
            return binarySearch(mid + 1 , high);

         return binarySearch(low , mid -1);
      }

      return 0;
   };


   return binarySearch( 0 , static_cast<int>(p.size()) - 1);


}

vector<double> Random :: UnitVersor(){

   vector<double> v(3);

   double theta = acos(1 - 2*Rannyu());
   double phi = Rannyu( 0 , 2*M_PI);

   v[0] = sin(theta)*cos(phi);
   v[1] = sin(theta)*sin(phi);
   v[2] = cos(theta);

   return v;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}



double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
