#ifndef __MATHOBJ__
#define __MATHOBJ__

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;
class MATDBL;
inline int mod(int i,int j){ return(i-(i/j)*j);}
inline double dmax1(double a,double b){return((a>b)?a:b);}
inline double dmin1(double a,double b){return((a>b)?b:a);}
//======================= VECINT =======================
class VECINT{
	public:
	int n;
	int* v;
	VECINT();
	VECINT(int);
	inline int& VECINT::operator[](int i){return v[i];}
	~VECINT();
};
//======================= VECINT =======================
class VECDBL{
	friend VECDBL operator*(MATDBL&,VECDBL&);
	public:
	int n;
	double* v;
	VECDBL();
	VECDBL(int);
	VECDBL(VECDBL&);
	void init(int);
	VECDBL(int,double);
	void init(int,double);
	VECDBL(int,double*);
	void init(int,double*);
    double norme(); 
	void operator=(const VECDBL&);
	void operator=(double c);
	VECDBL operator+(VECDBL&);
    VECDBL operator-(VECDBL&);
    double  operator*(VECDBL&);
    friend VECDBL operator*(double,VECDBL&);
    friend VECDBL operator*(VECDBL&,double);
	void operator+=(VECDBL&);
	void operator-=(VECDBL&);
    void operator*=(double);
    void operator/=(double);
	inline double& operator[](int i){return v[i];}
	~VECDBL();
};
//======================= MATINT =======================
class MATINT   {
public:
	int n,m;
	int **mat;
	MATINT();
	MATINT(int,int);
	void operator=(int c);
	inline int* operator[](int i){ return(mat[i]); }
	~MATINT();
};

class MATDBL   {
   friend  VECDBL operator*(MATDBL&,VECDBL&);
   public:
   int n ;
   int m; 
   double **mat;
	MATDBL();
   MATDBL(int,int);
	void init(int,int);
	MATDBL(MATDBL&);
   void operator=(const MATDBL&);
	friend void operator+=(const MATDBL&,const MATDBL& a);
	friend void operator-=(const MATDBL&,const MATDBL& a);
	void operator=(double c);
   void print();
	inline double* operator[](int i){ return(mat[i]); }
   friend MATDBL operator+(const MATDBL&,const MATDBL&); 
   friend MATDBL operator-(const MATDBL&,const MATDBL&); 
   friend MATDBL operator*(const MATDBL&,const MATDBL&);
   friend MATDBL operator+(const MATDBL&,double);
   friend MATDBL operator*(const MATDBL&,double);
   friend MATDBL operator*(double,const MATDBL&);
   MATDBL transpose();
   MATDBL inverse(double& determ);
   void solver_gauss(VECDBL&,VECDBL&);
   void solver_gaussLU(int,VECDBL&);
   ~MATDBL();
};
class TRIDBL   {
public:
    int n;  // nombre de lignes
    int m;  // nombre de colonnes`
    int k;  // 3ème dimension
    int s;  // totoal dimension
    double *tri;
	TRIDBL();
	TRIDBL(int,int,int);
	inline double& operator()(int a,int b, int c){return(tri[a*m*k+b*k+c]); }
	~TRIDBL();
};
 
#endif	
