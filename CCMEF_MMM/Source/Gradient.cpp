#include "Objet_math.h"
void SELEMbarre2d(VECDBL& vdle,MATDBL& vcore,VECDBL& vmate,MATDBL& sige);
void SELEMplanet3(VECDBL& vdle,MATDBL& vcore,VECDBL& vmate,MATDBL& sige);
//
//-----------------------------------------------------------------------
void gradient(int type,int nelt,int ndle,int ndln,int nnel,int ndim,int nsig,int nphy,
			  VECDBL& ddl,MATDBL& vmat,VECINT& kmat,MATDBL& coord,
			  MATINT& nelem,TRIDBL& sigma) {
//=======================================================================
//  calcul des contraintes aux noeuds selon le type d'éléments
//-----------------------------------------------------------------------
	int i,j,k,n,ieq,m,km,iel;
	VECDBL vdle(ndle),vmate(nphy);
	MATDBL sige(nnel,nsig),vcore(nnel,ndim);
	for (iel=0;iel<nelt;iel++) {
		km=kmat[iel]; for (i=0;i<nphy;i++) vmate[i]=vmat[km-1][i];
		m=0;
		for (j=0;j<nnel;j++) {
			n=nelem[iel][j];
			for(k=0;k<ndim;k++) vcore[j][k]=coord[n-1][k];
			for(k=0;k<ndln;k++) {
				ieq=(n-1)*ndln+k;
				vdle[m++]=ddl[ieq];
			}
		}
		switch(type) {
			case 1  : SELEMbarre2d(vdle,vcore,vmate,sige); break;
			case 2  : SELEMplanet3(vdle,vcore,vmate,sige); break;
			default : cout << "L'élément du type  "<<type<<"  n'est pas implanté dans gradient"<<endl;
		}
		for (j=0;j<nnel;j++) for(k=0;k<nsig;k++) sigma(iel,j,k)=sige[j][k];
	}
}

//-----------------------------------------------------------------------
void SELEMbarre2d(VECDBL& vdle,MATDBL& vcore,VECDBL& vmate,MATDBL& sige) 
{
	//=======================================================================
	//     calcul des contraintes axiales
	//     barre 2d 
	//-----------------------------------------------------------------------
	double young=vmate[0];
	double dx=vcore[1][0]-vcore[0][0];
	double dy=vcore[1][1]-vcore[0][1];
	double xlong=sqrt(dx*dx+dy*dy);
	double c=dx/xlong;
	double s=dy/xlong;
	double u1=vdle[0];
	double v1=vdle[1];
	double u2=vdle[2];
	double v2=vdle[3];
	sige[0][0]=sige[1][0]=young*((u2-u1)*c+(v2-v1)*s)/xlong;
}

//-----------------------------------------------------------------------
void SELEMplanet3(VECDBL& vdle,MATDBL& vcore,VECDBL& vmate,MATDBL& sige) 
{
//=======================================================================
//     calcul des contraintes aux noeuds
//     élasticité 2d (t3)
//-----------------------------------------------------------------------
}
