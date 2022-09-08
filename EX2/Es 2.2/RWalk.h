#ifndef __RWalk__
#define __RWalk__

#include <vector> 
#include <algorithm> 

#include "random.h"

using namespace std;

class RWalk
{
	public:

    //constructor
	RWalk(double a, int d , vector<double> O);


    //variables
    double get_dim();
	vector<double> get_coords();
	double get_step();
	double get_norm();

    //moves
	void move_to_origin();
	void RW_move_discrete();
	void RW_move_continuos(); //should be used only for dim=3 case, it is not general
		
	private:
	int dim;
	double step; 
	vector<double> x;
    vector<double> x0;
    Random rnd;
};

#endif