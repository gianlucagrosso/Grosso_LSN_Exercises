#include "Population.h"
#include <algorithm>
#include <cmath>
#include <iomanip>

void Population::NewIndividual (Individual ind){
	m_pop.push_back(ind);
}
long unsigned int Population::GetNIndividual(){ 
	return m_pop.size();
}
Individual Population::GetIndividual(int iInd){
	return m_pop[iInd];
}
bool Population::Check(){
	bool ris=true;
	for(long unsigned int iInd=0; iInd<m_pop.size(); iInd++) {
		ris = ris && m_pop[iInd].Check();
	}
	return ris;
}
void Population::StartingPop(){
	for(long unsigned int iInd=1; iInd<m_pop.size(); iInd++){
		m_pop[iInd].Shuffle();
	}
}
void Population::SortPop(){
	sort(m_pop.begin(),m_pop.end(),[](Individual a, Individual b){
		return a.CostFun() < b.CostFun();
          	});
}
void Population::PrintPop(const char* Filename){
	ofstream fout (Filename,ios::app);
	for(long unsigned int iInd=0; iInd<m_pop.size(); iInd++) {
        	//cout << "Individual " << iInd + 1 << endl;
		fout << "Individual " << iInd + 1 << endl;
                //m_pop[iInd].PrintInd();
                //cout << "La funzione costo L = "<< m_pop[iInd].CostFun() << endl << endl;
		fout << "La funzione costo L = "<< m_pop[iInd].CostFun() << endl << endl; 
        }
	fout.close();
}
void Population::PrintBest(int iInd,const char* Filename,int igen){
	ofstream fout (Filename,ios::app);
	int wd=12;
	fout << setw(wd) << igen << setw(wd) << m_pop[iInd].CostFun() << endl;
	fout.close();
}
void Population::Print_The_Best_Path(int iInd,const char* Filename){
	m_pop[iInd].PrintInd(Filename);
}
void Population::PrintHalfAverage(const char* Filename,int nindividual,int igen){
	ofstream fout (Filename,ios::app);
	int wd=12;
	int M = nindividual/2; //half of the population
	double sum=0., sum2=0., ave=0., ave2=0., err=0.; //variables
	for (int iInd=0; iInd<M; iInd++){
		sum += m_pop[iInd].CostFun();
		sum2 += pow(m_pop[iInd].CostFun(),2);
	}
	ave = sum/M; //averaged value
	ave2 = sum2/M; //square averaged value
	err = sqrt(fabs(ave2 - pow(ave,2))); //statistical uncertainty

	fout << setw(wd) << igen << setw(wd) << ave << setw(wd) << err << endl;
	fout.close();
}
int Population::Selection(int nind, Random* rnd, double exp){
	double r = rnd->Rannyu(); 
	int jInd = (int)(nind*pow(r,exp));
	return jInd; 
}
void Population::Mutation(int iInd, int whichmutation, int ncity,Random* rnd){
	if(whichmutation==1) m_pop[iInd].Gen_Mut_1((int)(ncity*(rnd->Rannyu())),(int)(ncity*(rnd->Rannyu())));
	if(whichmutation==2) m_pop[iInd].Gen_Mut_2((int)((ncity)*(rnd->Rannyu())),(int)((ncity)*(rnd->Rannyu())),(int)((ncity -1)*(rnd->Rannyu())));
	if(whichmutation==3) m_pop[iInd].Gen_Mut_3((int)(ncity*(rnd->Rannyu())),rnd,(int)((ncity/2)*(rnd->Rannyu())));
	if(whichmutation==4) m_pop[iInd].Gen_Mut_4((int)(ncity*(rnd->Rannyu())),(int)((ncity)*(rnd->Rannyu())));
}	
void Population::CrossOver(int cutindex,int Indmot,int Indfat){ // cut index must be in [1,Individual.size())
		//cout << "cutindex: " << cutindex << endl;
		Individual son1,son2;
		for(int icity=0; icity < cutindex; icity++) {
                	son1.NewCity(m_pop[Indmot].GetCity(icity));
                	son2.NewCity(m_pop[Indfat].GetCity(icity));
		}
		for(long unsigned int icity=0; icity<m_pop[Indmot].GetNcity(); icity++){
                	if(!son2.Found(0, cutindex, (m_pop[Indmot].GetCity(icity)).GetIndex())){
				son2.NewCity(m_pop[Indmot].GetCity(icity));
                	}
                	if(!son1.Found(0, cutindex, (m_pop[Indfat].GetCity(icity)).GetIndex())) {
				son1.NewCity(m_pop[Indfat].GetCity(icity));
                	}
		}
		m_pop[Indmot].DeleteFromIndex(cutindex);
		m_pop[Indfat].DeleteFromIndex(cutindex);
		for(long unsigned int icity=cutindex; icity < son1.GetNcity(); icity++) {
                	m_pop[Indmot].NewCity(son1.GetCity(icity));
                	m_pop[Indfat].NewCity(son2.GetCity(icity));
		}
}


