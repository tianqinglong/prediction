#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <R.h>
#include <Rinternals.h>

#include "overall_type2.h"
#include "brent.h"

int r=30, n=50;
double weiEta=1.5, weiBeta=1.5;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX];

SEXP transferDataToR(){
	int i;
	double betaMLE, etaMLE;
	double t, machep;

	SEXP MLEs_Data = PROTECT(allocVector(REALSXP,(4+n)));

	set_seed(time(NULL),533);
	simuDataType2(cenArray, realArray);

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.1,10,machep,t,deriEquation);

	etaMLE = 0;
	for(i=0;i<n;i++){
		etaMLE += pow(cenArray[i], betaMLE);
	}
	etaMLE = pow(etaMLE/r,1/betaMLE);

	REAL(MLEs_Data)[0] = n;
	REAL(MLEs_Data)[1] = r;

	REAL(MLEs_Data)[2] = betaMLE;
	REAL(MLEs_Data)[3] = etaMLE;

	for(i=0;i<n;i++){
		REAL(MLEs_Data)[i+4] = cenArray[i];
	}

	UNPROTECT(1);
	return MLEs_Data;
}
