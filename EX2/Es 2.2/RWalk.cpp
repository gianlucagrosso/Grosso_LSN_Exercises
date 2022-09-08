#include "RWalk.h"

RWalk::RWalk(double a, int d , vector<double> O)
{
	dim = d;
	step = a;
    x.resize(d);
    if( static_cast<int>(O.size()) == d ){
        x0 = O;
    } else cerr << "Origin vector dimension is different from RWalk dimension" << endl;
}

vector <double> RWalk::get_coords()
{
	return x;
}

double RWalk::get_step()
{
	return step;
}

void RWalk::move_to_origin()
{
	x = x0;
}

void RWalk::RW_move_discrete()
{
	vector<double> p ( 2*dim , static_cast<double>(1)/(2*dim));
    int dir=rnd.Discrete(p);
    	
	if(dir <= dim - 1)
        x[dir]+=step;
    else x[dir - dim]-=step;
}

void RWalk::RW_move_continuos(){
    if (dim == 3){

        vector<double> v = rnd.UnitVersor();

        for(int i = 0; i < dim; i++)
            x[i] += step*v[i];


    }else cerr << "Random Walk dimension is >3" << endl;


}


double RWalk::get_norm()
{
	double norm=0;
	
	for (int i = 0; i < dim; i++)
		norm += x[i]*x[i] - x0[i]*x0[i];

	return norm;
}