#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <R.h>
#include <Rinternals.h>

#include "brent.h"
#define ARRAY_MAX 100

int r,n;
double Pt=0.25, Er=10;
double weiEta=1.5, weiBeta=1.5;

double cenArray[ARRAY_MAX], comArray[ARRAY_MAX];

double deriEquation1(double weibullBeta);
void simuDataType1(double *cenVec, double *comVec);

SEXP transToR1(){
	int i;
	double t, machep;
	double betaMLE, etaMLE;

	set_seed(time(NULL),500);
	r=0;
	while(r < 2){
		simuDataType1(cenArray, comArray);
	}

	SEXP MLEs_Data = PROTECT(allocVector(REALSXP,(4+n)));

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.1,10,machep,t,deriEquation1);
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