#include "Objet_math.h"
void modifKF(int methode,MATINT& indice, MATDBL& valimp,MATDBL& vkg,VECDBL& vfg) {
//======================================================================
// imposition des conditions aux limites
//======================================================================
	int i, j, k, ieq;
	int nnt = indice.n, ndln = indice.m, neq=vfg.n;
	if (methode == 1) {
		// utilisation de la methode du terme diagonal dominant
		double alfa = 1.0e40;
		for (i = 0; i < nnt; i++) {
			for (j = 0; j < ndln; j++) {
				if (indice[i][j] == 1) {
					ieq = i*ndln + j;
					vkg[ieq][ieq] += alfa;
					vfg[ieq] = valimp[i][j] * alfa;
				}
			}
		}
	}
	else if (methode==2) {
		// utilisation de la methode du terme diagonal unité
		for (i = 0; i<nnt; i++) {
			for (j = 0; j<ndln; j++) {
				if (indice[i][j] == 1) {
					ieq = i*ndln + j;
					for (k = 0; k < neq; k++) vfg[k] -= vkg[ieq][k] * valimp[i][j];
					for (k = 0; k < neq; k++) vkg[k][ieq] = vkg[ieq][k] = 0;
					vkg[ieq][ieq] = 1.0;
					vfg[ieq] = valimp[i][j];
				}
			}
		}
	}
}