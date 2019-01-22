#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall_type2.h"

void simuDataType2(double* cenVec, double* comVec){
	int i;
	double curOrd, u, nextOrd;

	curOrd=0;

	for(i=0;i<n;i++){
		u = unif_rand();
		nextOrd = 1 - (1-curOrd)*pow(1-u, 1.0/(n-i));
		comVec[i] = weiEta*pow(-log(1-nextOrd), 1/weiBeta);
		cenVec[i] = (i<r)?comVec[i]:comVec[r-1];
		curOrd = nextOrd;
	}
}