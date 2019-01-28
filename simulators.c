#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall.h"

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

void simuDataType2(double* cenVec, double* comVec, double shape, double scale){
    int i;
    double curOrd, u, nextOrd;
    
    curOrd=0;
    
    for(i=0;i<n;i++){
        u = unif_rand();
        nextOrd = 1 - (1-curOrd)*pow(1-u, 1.0/(n-i));
        comVec[i] = scale*pow(-log(1-nextOrd), 1/shape);
        cenVec[i] = (i<r)?comVec[i]:comVec[r-1];
        curOrd = nextOrd;
    }
}
