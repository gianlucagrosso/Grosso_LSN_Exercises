#include "Individual.h"
#include "funzioni.h" //Pbc
#include <algorithm>
#include <iomanip>

void Individual::NewCity (double x, double y, int j){
	City newcity(x,y,j);
	m_Ind.push_back(newcity);
}
void Individual::NewCity (const City & c){
            City temp(c.GetX(), c.GetY(), c.GetIndex());
            m_Ind.push_back(temp);
}
City Individual::GetCity(int icity){
	return m_Ind[icity];
}
void Individual::SetCity (int i, double x, double y, int j ){
	m_Ind[i].SetX(x);
	m_Ind[i].SetY(y);
	m_Ind[i].SetIndex(j);
}
void Individual::SetCity (int i, const City & c) {
            m_Ind[i].SetX(c.GetX());
            m_Ind[i].SetY(c.GetY());
            m_Ind[i].SetIndex(c.GetIndex());
        }

bool Individual::Check(){
	bool ris=true;
        vector<int> labels;
	for(long unsigned int icity=0; icity<m_Ind.size(); icity++){
		labels.push_back(m_Ind[icity].GetIndex());
	}
	ris = ris && (m_Ind[0].GetIndex()==1);
	for(long unsigned int icity=1; icity<=m_Ind.size(); icity++){
		ris = ris&&(find(labels.begin(), labels.end(), icity) != labels.end());
	}
	return ris;
}
void Individual::Shuffle(){
	random_shuffle(m_Ind.begin()+1,m_Ind.end());
}
double Individual::CostFun(){
	double L = 0;
	int ncity = (int) m_Ind.size();
	for (long unsigned int icity=0; icity < m_Ind.size(); icity++){
		L += m_Ind[icity].Distanza(m_Ind[Pbc(icity+1,ncity)]);
	}
	return L;
}
bool Individual::Found(int index1, int index2, int label){
	vector<int> labels;
	for(int icity=index1; icity<index2; icity++){
		labels.push_back(m_Ind[icity].GetIndex());
	}
	return (find(labels.begin(), labels.end(), label) != labels.end());
}
void Individual::Gen_Mut_1(int i,int j){ 
	City City1;
	City1 = m_Ind[0];
	m_Ind.erase(m_Ind.begin()); // remove the first in order to preserve the position of the first city
	int ncity = (int) m_Ind.size();  // ncity - the first                       
	iter_swap(m_Ind.begin()+ Pbc(i,ncity), m_Ind.begin()+ Pbc(j,ncity)); 
	m_Ind.insert(m_Ind.begin(), City1);  //after mutation reinsert the first city in the position 1 of the path
}

void Individual::Gen_Mut_2(int index, int shift,int contigcity){ //contigcity must be in [0,(m_Ind.size())-1)
		
	City City1;
	City1 = m_Ind[0];
	m_Ind.erase(m_Ind.begin()); // remove the first in order to preserve the position of the first city
	int ncity = (int) m_Ind.size();  // ncity - the first  

	for(int i = 0; i< shift; i++){	//repeat the move shift times
		for(int j = 0 ; j < contigcity; j++ ){ //move the contigcity back one place
			iter_swap(m_Ind.begin() + Pbc(index-i+j+ncity, ncity) , m_Ind.begin() + Pbc(index-i+j+ncity-1, ncity) ); //+ncity is to deal with negative numbers

		}
	}
	m_Ind.insert(m_Ind.begin(), City1);  //after mutation reinsert the first city in the position 1 of the path	

}

void Individual::Gen_Mut_3(int index1 ,Random* rnd, int contigcity){  //contigcity must be in [0,(m_Ind.size())/2-1]
	City City1;
	City1 = m_Ind[0];
	m_Ind.erase(m_Ind.begin()); // remove the first in order to preserve the position of the first city
	int ncity = (int) m_Ind.size(); // ncity - the first
	int index2 = (int)((m_Ind.size()-2.*contigcity + 1)*(rnd->Rannyu())); // to avoid ovrelaps between the m contiguos city

	for (int icontig = 0; icontig< contigcity; icontig++){
		iter_swap (m_Ind.begin()+ Pbc(index1 + icontig,ncity),m_Ind.begin()+ Pbc(index1 + contigcity + index2 + icontig,ncity));
	}
	
	m_Ind.insert(m_Ind.begin(), City1);  //after mutation reinsert the first city in the position 1 of the path	
}

void Individual::Gen_Mut_4(int index ,int contigcity){ //contigcity must be in  [0,m_ind.size()-1]
	City City1;
	City1 = m_Ind[0];
	m_Ind.erase(m_Ind.begin());// remove the first in order to preserve the position of the first city
	int ncity = (int) m_Ind.size(); // ncity - the first

	if ((long unsigned int)(index + contigcity) > m_Ind.size()){
		rotate(m_Ind.begin(),m_Ind.begin() + Pbc(index ,ncity) ,m_Ind.end());
		reverse( m_Ind.begin(), m_Ind.begin() + contigcity );
		rotate(m_Ind.begin(),m_Ind.end() - index ,m_Ind.end());
	}
	else{
		reverse( m_Ind.begin() + index, m_Ind.begin() + index + contigcity);
	}
	m_Ind.insert(m_Ind.begin(), City1); //after mutation reinsert the first city in the position 1 of the path	
}
void Individual::PrintInd(const char* Filename){
	ofstream fout (Filename,ios::app);
	int wd=12;
	for(long unsigned int icity=0; icity < m_Ind.size(); icity++) {
        	fout << setw(wd) << m_Ind[icity].GetIndex() << setw(wd) << m_Ind[icity].GetX() << setw(wd) << m_Ind[icity].GetY() << endl;
	}
	fout.close();
}



	
