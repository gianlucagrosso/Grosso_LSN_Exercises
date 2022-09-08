#include "funzioni.h"

int Pbc(int icity,int ncity){
	//int index;
	//cout << "L'indice inserito : " << icity << endl;
	if (icity >= ncity){
		//cout << "L'indice modificato : " << icity-ncity << endl << endl;
		return (icity%ncity);
	}
	else{
		//cout << "L'indice modificato : " << icity << endl << endl;
		return icity;
	}
}
