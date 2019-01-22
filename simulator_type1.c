#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall_type1.h"

void simuDataType1(double *cenVec, double *comVec){
	int i;
	double curOrd, nextOrd, u;

	n = Er/Pt;

	r = 0;
	curOrd = 0;
	for(i=0; i<n; i++){
		u = unif_rand();
		nextOrd = 1 - (1-curOrd)*pow(1-u,1.0/(n-i));
		r += (nextOrd <= Pt);
		comVec[i] = weiEta*pow(-log(1-nextOrd), 1/weiBeta);
		cenVec[i] = comVec[i];
		curOrd = nextOrd;
	}
	for(i=r;i<n;i++){
		cenVec[i] = weiEta*pow(-log(1-Pt), 1/weiBeta);
	}
}