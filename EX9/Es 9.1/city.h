#ifndef __city_h__
#define __city_h__

#include <iostream>

using namespace std;

class City {
	public:	
		City(); //default constructor
		City(double,double,int); //constructor 
		
		City(const City&); //copy constructor
		~City(){;}; //distructor

		//cartesian coordinate
		double GetX() const;
		double GetY() const;
		int GetIndex() const;

		void SetX(double);
		void SetY(double);
		void SetIndex(int);

		//polar coordinate
		double GetR() const;
		double GetTheta() const;

		double Distanza(const City&) const; //distance from point 

		bool operator == (const City&) const;
		
	private:
		double m_x, m_y; //plane coordinate
		int m_index; //index that labels each city generated
};
#endif
