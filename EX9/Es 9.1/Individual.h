#ifndef __Individual_h__
#define __Individual_h__

#include <fstream>
#include <vector>
#include "city.h"
#include "random.h"

using namespace std;


class Individual{
		
		public:
			~Individual(){;}; //distructor

			void NewCity (double,double,int); //constructor adding a new city
			void NewCity (const City &);

			long unsigned int GetNcity() { return m_Ind.size();};
			City GetCity(int);

			void SetCity (int,double,double,int); //modify the i-th city
			void SetCity (int,const City &);

			void DeleteFromIndex(int index) {m_Ind.erase(m_Ind.begin() + index, m_Ind.end());};
			void Clear(){m_Ind.clear();}; //clear the entire vector

			bool Check(); //check function which verify if every individuals fulfils the bonds

			void Shuffle(); //shuffle randomly the city in Individuals

			double CostFun(); //return the length og the path: the Cost Function

			bool Found(int,int,int); //verify if finds city's label in the Individual

			//genetic mutation
			void Gen_Mut_1(int,int);
			void Gen_Mut_2(int, int, int);
			void Gen_Mut_3(int,Random*,int);
			void Gen_Mut_4(int,int);

			void PrintInd(const char* Filename); //print the Individual			

		private:
			vector <City> m_Ind; //Individual: cities' contenitor

};
#endif
