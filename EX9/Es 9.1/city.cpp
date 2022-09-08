#include "city.h"
#include <cmath>

City::City(){
	m_x = 0.;
	m_y = 0.;
	m_index = 0;
}
City::City(double x, double y, int i){
	m_x = x;
	m_y = y;
	m_index = i;
}
City::City(const City& c){
	m_x = c.GetX();
	m_y = c.GetY();
	m_index = c.GetIndex();	
}
double City::GetX() const{
	return m_x;
}
double City::GetY() const{
	return m_y;
}
int City::GetIndex() const{
	return m_index;
}
void City::SetX(double x){
	m_x = x;
}
void City::SetY(double y){
	m_y = y;
}
void City::SetIndex(int i){
	m_index = i;
}
double City::GetR() const{
	return sqrt(pow(m_x,2) + pow(m_y,2));
}
double City::GetTheta() const{
	return atan(GetY()/GetX());
}
double City::Distanza(const City& c) const{
	return sqrt (pow(GetX()-c.GetX(),2) + pow(GetY()-c.GetY(),2));
}
bool City::operator==(const City &c) const {
    return (m_x == c.GetX()) && (m_y == c.GetY());
}
