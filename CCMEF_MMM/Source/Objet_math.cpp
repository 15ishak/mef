#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Objet_math.h"
using namespace std;
/******************** VECINT *******************/
VECINT::VECINT(){n=0;v=NULL;}
VECINT::VECINT(int a){
	n=a;v=new int[n];
	for(int i=0;i< n;i++) v[i]=0;
}
VECINT::~VECINT(){ if(n!=0) delete v;}
/******************** VECDBL *******************/
VECDBL::VECDBL(){n=0;v=NULL;}
VECDBL::VECDBL(int a){n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=0.;}
VECDBL::VECDBL(VECDBL& vv){
   n=vv.n;v=new double[n];
   double *pvv,*pv;
   pv=v;
   pvv=vv.v;for(int i=0;i< n;i++)*pv++=*pvv++;
}
void VECDBL::init(int a){
   delete[]  v;n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=0.;
}
VECDBL::VECDBL(int a,double c){
   n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=c;
}
void VECDBL::init(int a,double c){
   delete[]  v;n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=c;
}
VECDBL::VECDBL(int a,double* c){
   n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=c[i];
}
void VECDBL::init(int a,double* c){
   delete[]  v;n=a;v=new double[n];for(int i=0;i< n;i++)v[i]=c[i];
}
double VECDBL::norme(){
   double xn=0.;int i;
   for(i=0;i< n;i++)xn+=v[i]*v[i];
   xn=sqrt(xn);
   return(xn);
}
void VECDBL::operator=(const VECDBL& vv){
   double *pvv,*pv; 
   pv=v;
   pvv=vv.v;
   if(vv.n == n) for(int i=0;i< n;i++) *pv++=*pvv++;
   else
   cout<<"*** Erreur : operateur = entre deux vecteurs de tailles differentes\n";
}
void VECDBL::operator=(double c){for(int i=0;i< n;i++)v[i]=c;}
VECDBL VECDBL::operator+(VECDBL& vv){
   VECDBL s(n);
   double *pvv,*pv,*ps;
   pv=v;pvv=vv.v;ps=s.v;
   if(vv.n == n){for(int i=0;i< n;i++){*ps++ = *pvv++ + *pv++;}return(s);}
   else
   {cout<<"*** Erreur : operateur + entre deux vecteurs de tailles differentes\n";exit(0);}
	return(s);
}
VECDBL VECDBL::operator-(VECDBL& vv){
   VECDBL s(n); 
   double *pvv,*pv,*ps; 
   pv=v;pvv=vv.v;ps=s.v;
   if(vv.n == n){for(int i=0;i< n;i++){*ps++=*pv++ - *pvv++;}return(s);} 
   else 
   {cout<<"*** Erreur : operateur - entre deux vecteurs de tailles differentes\n";exit(0); }
	return(s);
}
VECDBL operator*(double d, VECDBL& vec){
   VECDBL s(vec.n);
   double *pv,*ps;
   pv=vec.v;ps=s.v;
   for(int i=0;i< vec.n;i++){*ps++ = *pv++ * d;}return(s);
}
VECDBL operator*(VECDBL& vec,double d){VECDBL x(vec.n); x=d*vec; return(x);}
double VECDBL::operator*(VECDBL& vv){
   double d=0.0;
   double *pvv,*pv;
   pv=v;pvv=vv.v;
   if(vv.n == n){for(int i=0;i< n;i++){d += *pvv++ * *pv++;} return(d);}
   else    
   {cout<<"*** Erreur : produit scalaire entre deux vecteurs de tailles differentes\n";exit(0);}
	return(d);
}
VECDBL operator*(MATDBL& m, VECDBL& v){
   VECDBL s(m.n);
   if(v.n == m.m )
   {for(int i=0;i < m.n;i++){ s.v[i]=0. ;
         for(int j=0;j<m.m;j++)s.v[i]+=m.mat[i][j]*v.v[j];
      }
      return(s);
   }
   else {cout<<"*** Erreur : produit matrice par vecteur de tailles incompatibles";exit(0);}
   return(s);
}
void VECDBL::operator+=(VECDBL& vv){
   double *pvv,*pv;
   pv=v;pvv=vv.v;;
   if(vv.n == n)for(int i=0;i< n;i++){*pv++ += *pvv++ ;}
   else    
   {cout<<"*** Erreur : operateur += entre deux vecteurs de tailles differentes\n";}
}
void VECDBL::operator-=(VECDBL& vv){
   double *pvv,*pv;
   pv=v;pvv=vv.v;;
   if(vv.n == n)for(int i=0;i< n;i++){*pv++ -= *pvv++ ;}
   else
   {cout<<"*** Erreur : operateur -= entre deux vecteurs de tailles differentes\n";}
}
void VECDBL::operator*=(double d){
   double *pv;   pv=v;  for(int i=0;i< n;i++){*pv++ *= d ;}
}
void VECDBL::operator/=(double d){ 
   double *pv;   pv=v;  for(int i=0;i< n;i++){*pv++ /= d ;}
}
VECDBL::~VECDBL(){ delete[] v;}
/******************** MATINT *******************/
MATINT::MATINT(){n=0;m=0;mat=NULL;}
MATINT::MATINT(int a, int b) {
	n=a;m=b;
	typedef int* pdb;
	mat=new pdb[n];
	for(int i=0;i<n;i++) {
		mat[i]=new int[m];
		if(!mat[i]) {cout<< "*** Erreur : plus d'esapace memoire"<<"\n";exit(0);}
	}
	for(int i=0;i<n;i++) for(int j=0;j<m;j++) mat[i][j]=0;
}
void MATINT::operator=(int c){
	int i,j;
	for(i=0;i<n;i++) for(j=0;j<m;j++) mat[i][j]=c;
}
MATINT::~MATINT() {
	if(n!=0 || m!=0){
		for(int i=0;i<n;i++) delete[] mat[i];
		delete mat;
	}
}
/***********************************************************************/
/***********************************************************************/
MATDBL::MATDBL(){n=0;m=0;mat=0;}
MATDBL::MATDBL(int a, int b)
{
   n=a;m=b;
   typedef double* pdb;
   mat=new pdb[n];
   for(int i=0;i<n;i++) {mat[i]=new double[m];
      if(!mat[i]) cout<< "plus d'esapace memoire"<<"\n";
   }
   for(int i=0;i<n;i++) for(int j=0;j<m;j++) mat[i][j]=0.0;
}
void MATDBL::init(int a, int b)
{
	int i,j;
   for(i=0;i<n;i++) delete[] mat[i];
   delete[] mat;
   n=a;m=b;
   typedef double* pdb;
   mat=new pdb[n];
   for(i=0;i<n;i++) mat[i]=new double[m];
   for(i=0;i<n;i++) for(j=0;j<m;j++) mat[i][j]=0.0;
}
MATDBL::MATDBL(MATDBL& a){
   n=a.n;m=a.m;
   typedef double* pdb;
   mat = new pdb[n];
   for(int i=0;i<n;i++) mat[i]=new double[m];
   for(int i=0;i<n;i++) for(int j=0;j<m;j++) mat[i][j]=a.mat[i][j];
}
void MATDBL::print()
{
   cout<<"\n";
   for(int i=0;i<n;i++)
   {
      for(int j=0;j<m;j++)
      { cout<<mat[i][j]<<"\t";}
      cout<<endl;
   }
}
void operator+=(const MATDBL& a,const MATDBL& b){
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++) a.mat[i][j]+=b.mat[i][j];
}
void operator-=(const MATDBL& a,const MATDBL& b){
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++) a.mat[i][j]-=b.mat[i][j];
}
void MATDBL::operator=(const MATDBL& b)
{
   for(int i=0;i<n;i++) for(int j=0;j<m;j++) mat[i][j]=b.mat[i][j];
}
void MATDBL::operator=(double c){
   for(int i=0;i<n;i++) for(int j=0;j<m;j++) mat[i][j]=c;
}
MATDBL operator*(const MATDBL& a,double d){
   MATDBL r(a.n,a.m);
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++)r.mat[i][j]= a.mat[i][j]*d;
   return(r);
}
MATDBL operator*(double d,const MATDBL& a){
   MATDBL r(a.n,a.m);
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++)r.mat[i][j]= a.mat[i][j]*d;
   return(r);
}
MATDBL operator+(const MATDBL& a,const MATDBL& b)
{
   MATDBL r(a.n,a.m);
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++)
   r.mat[i][j]=a.mat[i][j]+b.mat[i][j];
   return(r);
}
MATDBL operator-(const MATDBL& a,const MATDBL& b)
{
   MATDBL r(a.n,a.m);
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++)
   r.mat[i][j]=a.mat[i][j]-b.mat[i][j];
   return(r);
}
MATDBL operator*(const MATDBL& a,const MATDBL& b)
{
   MATDBL r(a.n,b.m);
   if( a.m != b.n){
      cout<<" calcul de produit de deux matrices "<< a.n<<"X"<<a.m;
      cout<< " et "<<b.n<<"X"<<b.m<<"  est impossible"<<"\n";
      exit(0);
   }
   double val;
   for(int i=0;i<a.n;i++)
   {
      for(int j=0;j<b.m ;j++)
      {
         val=0.;
         for(int k=0;k<a.m;k++) val+=a.mat[i][k]*b.mat[k][j];
         r.mat[i][j]=val;
      }
   }
   return(r);
}
MATDBL operator+(const MATDBL& a,double d)
{
   MATDBL r(a.n,a.m);
   for(int i=0;i<a.n;i++) for(int j=0;j<a.m;j++)
   r.mat[i][j]=a.mat[i][j]+d;
   return(r);
}
MATDBL MATDBL::transpose(){
   MATDBL b(m,n);
   for(int i=0;i<n;i++) for(int j=0;j<m;j++)
   b.mat[j][i] = mat[i][j];
   return(b);
}
MATDBL MATDBL::inverse(double& determ)
{
   MATDBL b(n,m);
   int i,j,k,l,s,nb_per;
   int* jj; jj=new int[n];
   double epsi=1.e-20;double ss;
   for(i=0;i<n;i++) for(j=0;j<n;j++) b.mat[i][j] = mat[i][j];
   double dd=0;
   double pd=1.;
   for(l=0;l<n;l++){
      dd=0.;
      for(k=0;k<n;k++)dd+=b.mat[l][k]*b.mat[l][k];
      dd=sqrt(dd); pd*=dd;
   }
   determ=1.;
   for(l=0;l<n;l++)jj[l]=l;
   double cc;
   nb_per=0;
   for(l=0;l<n;l++){
      cc=0.;s=l;
      for(k=l;k<n;k++){
         if( fabs(b.mat[l][k]) > fabs(cc)){
            s=k;cc=b.mat[l][k];
         }
      }
      if( l != s){
         nb_per++;
         k=jj[s];jj[s]=jj[l];jj[l]=k;
         for(k=0;k<n;k++){
            ss=b.mat[k][l];
            b.mat[k][l]=b.mat[k][s];
            b.mat[k][s]=ss;
         }
      }
      b.mat[l][l]=1.;
      determ=determ*cc;
      for(k=0;k<n;k++)b.mat[l][k]/=cc;
      for(k=0;k<n;k++){
         if( l != k){
            cc=b.mat[k][l];
            if(fabs(cc) >epsi){
               b.mat[k][l]=0;
               for(j=0;j<n;j++)b.mat[k][j]-=cc*b.mat[l][j];
            }
         }
      }
   }
   for(l=0;l<n;l++){
      if(jj[l] != l){
         s=l;s++;
         while ( s < n-1 && jj[s] != l)s++;
         jj[s]=jj[l];
         for(k=0;k<n;k++){
            cc=b.mat[l][k];
            b.mat[l][k]=b.mat[s][k];
            b.mat[s][k]=cc;
         }
         jj[l]=l;
      }
   }
   if(mod(nb_per,2)){determ=-determ;}
   delete[]  jj;
   return(b);
}
//
MATDBL::~MATDBL() {
   for(int i=0;i<n;i++) delete[] mat[i];
   delete[] mat;
}
void MATDBL::solver_gauss(VECDBL& b,VECDBL& x) {
   int ii,ij,ieq1,ij1;
   int neq=b.n;
   double zero=0.0;
   double cf;
   int nm1=neq-1;
   for (int ieq=0; ieq<nm1; ieq++) {
	  if (mat[ieq][ieq]==zero) {
		 cout<<"*** pivot nul **** equation no: "<<ieq+1<<"\n"<<flush;
		 exit(1);
	  }
	  ieq1=ieq+1;
	  for (ii=ieq1; ii<neq; ii++){
		cf=mat[ii][ieq];
		if (cf!=zero) {
		   cf=cf/mat[ieq][ieq];
		   b[ii]-=(cf*b[ieq]);
		   for (ij=ieq1; ij<neq; ij++) mat[ii][ij]-=(cf*mat[ieq][ij]);
		}
	  }
   }
   if (mat[nm1][nm1]==zero) {
	  cout<<"*** pivot nul **** equation no: "<<neq<<"\n"<<flush;
	  exit(1);
   }
   else x[nm1]=b[nm1]/mat[nm1][nm1];
   for (ii=0; ii<nm1; ii++) {
	  ieq1--;
	  cf=zero;
	  ij1=ieq1+1;
	  for (ij=ij1; ij<neq; ij++) cf+=(mat[ieq1][ij]*x[ij]);
	  x[ieq1]=(b[ieq1]-cf)/mat[ieq1][ieq1];
   }
}
void MATDBL::solver_gaussLU(int option,VECDBL& b) {
   int i,j,k,j1,kp1,np1;
   int neq=b.n;
   double zero=0.0;
   double cf;
   int nm1=neq-1;

	if (option==0 || option==1) {
	   for (int k=0; k<nm1; k++) {
		 if (mat[k][k]==zero) {
			cout<<"*** pivot nul **** equation no: "<<k+1<<"\n"<<flush;
			exit(1);
		 }
		 kp1=k+1;
		 for (i=kp1; i<neq; i++){
			cf=-mat[i][k]/mat[k][k];
			mat[i][k]=cf;
			for (j=kp1; j<neq; j++) mat[i][j]+=(cf*mat[k][j]);
		 }
	   }
   }
   if (option==0 || option==2) {
		np1=neq+1;

		for (k=0; k<nm1; k++) {
		  kp1=k+1;
		  for (i=kp1; i<neq; i++) b[i]+=(mat[i][k]*b[k]);
		}
		if (mat[neq-1][neq-1]<1.0e-20) {
		 cout<<"*** pivot nul **** equation no: "<<neq<<"\n"<<flush;
		 exit(1);
		}
		b[neq-1]/=mat[neq-1][neq-1];
		for (k=2; k<=neq; k++) {
		  i=np1-k-1;
		  if (mat[i][i]<1.0e-20) {
			 cout<<"*** pivot nul **** equation no: "<<k<<"\n"<<flush;
			 exit(1);
		  }
		  j1=i+1;
		  for (j=j1; j<neq; j++)  b[i]-=(b[j]*mat[i][j]);
		  b[i]/=mat[i][i];
		}
   }
}
/******************** TRIDBL *******************/
TRIDBL::TRIDBL(){n=0;m=0;k=0;s=0;tri=NULL;}
TRIDBL::TRIDBL(int a, int b, int c) 
{
	int i;
	n=a; m=b; k=c; s = n*m*k;
	tri=new double[s];
	if(!tri) {
		cout<< "plus d'esapace memoire"<<"\n";
		exit(0);
	}

	for(i=0;i<s;i++) tri[i]=(double)0.0;
}

TRIDBL::~TRIDBL() {
	if(s!=0) delete tri;
}
