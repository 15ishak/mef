#include <iomanip>
#include "Objet_math.h"
using namespace std;

void ELEMbarre2d(MATDBL& vcore,VECDBL& vmate,MATDBL& vke,VECDBL& vfe);
void ELEMplanet3(MATDBL& vcore,VECDBL& vmate,MATDBL& vke,VECDBL& vfe);
void elemlib(int type,MATDBL& vcore,VECDBL& vmate,MATDBL& vke,VECDBL& vfe) {
	switch(type) {
		case 1 : ELEMbarre2d(vcore,vmate,vke,vfe); break;
		case 2 : ELEMplanet3(vcore,vmate,vke,vfe); break;
		default : cout << "L'element type  "<<type<<"  n'est pas encore implante"<<endl; exit(0);
	}
}

void ELEMbarre2d(MATDBL& vcore,VECDBL& vmate,MATDBL& vke,VECDBL& vfe) {
//=======================================================================
//   calcul de la matrice de rigidité et du vecteur élémentaires
//   élément de barre 2D
//=======================================================================
	double dx=vcore[1][0]-vcore[0][0];
	double dy=vcore[1][1]-vcore[0][1];
	double xl=sqrt(dx*dx+dy*dy);
	if (xl==0.0) {
		cout<<" *** ELEMbarre2d : probleme des coordonnees dans l element "<<endl;
		exit(1);
	}
	double young=vmate[0];
	double sect=vmate[1];
	double rho=vmate[2];
	double vk=sect*young/xl;
	double c=dx/xl;
	double s=dy/xl;
	double c2=vk*c*c;
	double s2=vk*s*s;
	double cs=vk*c*s;
	vke[0][0] = vke[2][2] = c2;
	vke[0][2] = vke[2][0] = -c2;
	vke[1][1] = vke[3][3] = s2;
	vke[1][3] = vke[3][1] = -s2;
	vke[0][1] = vke[1][0] = vke[2][3] = vke[3][2] = cs;
	vke[1][2] = vke[2][1] = vke[0][3] = vke[3][0] = -cs;
	vfe[1] = vfe[3] = -4.9*rho*xl*sect;
}

void ELEMplanet3(MATDBL& vcore,VECDBL& vmate,MATDBL& vke,VECDBL& vfe) 
{
	//=======================================================================
	//   calcul de la matrice de rigidité et du vecteur élémentaires
	//   élément 2D - T3 
	//=======================================================================

}
