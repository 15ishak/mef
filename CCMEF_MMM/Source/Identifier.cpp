#include <iomanip>
#include <fstream>
using namespace std;
void identifier(int type,fstream& Out,int& idasc,char* dofname,
                int& ndim,int& ndln,int& nnel,int& nsig) {
	switch(type) {
		case 1  : 
			Out<<"\tType d'éléments : barre 2D"<<endl;
			idasc=22; ndim=2; ndln=2; nnel=2; nsig=1;
			strcpy_s(dofname,100,"5  UX   UY  Rx  Ry  sigm");
			break;
		case 2  : 
			Out<<"\tType d'éléments : plane T3"<<endl;
			idasc=23; ndim=2; ndln=2; nnel=3; nsig=4;
			strcpy_s(dofname,100,"8  UX   UY   Rx  Ry  Sxx  Syy  Sxy Mises");
			break;
		default : Out<<"L'element type  "<<type<<"  n'est pas encore implante"<<endl; exit(0);
	}
}