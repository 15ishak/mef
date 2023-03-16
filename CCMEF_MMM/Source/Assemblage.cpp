#include "Objet_math.h"

void assemblage(VECINT& loce, MATDBL& vke, VECDBL& vfe, MATDBL& vkg, VECDBL& vfg) {
	/*c=====================================================================
	c
	c     assemblage d'une MATDBL et/ou d'un VECDBL elementaire
	c       entrees
	c          ndle    nombre de d.l. de l'element
	c          loce    VECDBL de localisation de l'element
	c          vke     MATDBL elementaire ke(pleine)
	c          vfe     VECDBL elementaire fe
	c       sorties
	c          vkg     MATDBL globale (MATDBL pleine)
	c          vfg     VECDBL sollicitations global
	c=========================== debut du code executable ==================
	c
	c*/
	int i, j, ig, jg;
	int ndle = loce.n;
	for (i = 0; i<ndle; i++) {
		ig = loce[i];
		for (j = 0; j<ndle; j++) {
			jg = loce[j];
			vkg[ig][jg] += vke[i][j];
		}
		vfg[ig] += vfe[i];
	}
}
