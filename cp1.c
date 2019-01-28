#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <R.h>
#include <Rinternals.h>

#include "brent.h"
#include "overall.h"

int r, n;
double Er, Pt;
double weiEta, weiBeta, alpha;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX], weightArray[ARRAY_MAX];
double cenBootSample[ARRAY_MAX], comBootSample[ARRAY_MAX];
double etaMLE, betaMLE;
double betaBoot[B], etaBoot[B];

SEXP transToR1(SEXP Rseed, SEXP Rlower, SEXP Rupper, SEXP REr, SEXP RPt, SEXP Reta, SEXP Rbeta, SEXP Rboot){
	int i, j, seed = asInteger(Rseed), bootType = asInteger(Rboot);
	double t, machep, nt, dt;
	double betaMLEBoot, etaMLEBoot;
	double interGPQ, interPerBoot, interFonseca, interPlugIn;
	double lAlpha = asReal(Rlower), uAlpha = asReal(Rupper);

	Er = asReal(REr); Pt = asReal(RPt);
	weiEta = asReal(Reta); weiBeta = asReal(Rbeta);

	set_seed(time(NULL)+seed,seed);

	r=0;
	while(r < 2){
		simuDataType1(cenArray, realArray);
	}

	SEXP MLEs_Data = PROTECT(allocVector(REALSXP,(4+n+2*B+8)));

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.01,100,machep,t,deriEquation);
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

	if(bootType){
		for(j=0;j<B;j++){
			getWeight(weightArray);
			betaMLEBoot = zero(0.01,100,machep,t,deriFRWB);
			nt=0;dt=0;
			for(i=0;i<n;i++){
				nt += weightArray[i]*pow(cenArray[i],betaMLEBoot);
				dt += ((i<r) ? weightArray[i] : 0);
			}
			etaMLEBoot = pow(nt/dt, 1/betaMLEBoot);

			betaBoot[j] = betaMLEBoot;
			etaBoot[j] = etaMLEBoot;

			REAL(MLEs_Data)[4+n+j] = betaMLEBoot;
			REAL(MLEs_Data)[4+n+B+j] = etaMLEBoot;
		}
	}
	else{

	}

	alpha = lAlpha;
	interGPQ = zero(0.005,100,machep,t,getIntervalsGPQ);
	interPerBoot = zero(0.005,100,machep,t,getIntervalsPercentileBoot);
	interFonseca = zero(0.005,100,machep,t,getIntervalsFonseca);
	interPlugIn = qweibull(alpha, betaMLE, etaMLE, 1, 0);
	REAL(MLEs_Data)[4+n+2*B] = interGPQ;
	REAL(MLEs_Data)[4+n+2*B+1] = interPerBoot;
	REAL(MLEs_Data)[4+n+2*B+2] = interFonseca;
	REAL(MLEs_Data)[4+n+2*B+3] = interPlugIn;

	alpha = uAlpha;
	interGPQ = zero(0.005,100,machep,t,getIntervalsGPQ);
	interPerBoot = zero(0.005,100,machep,t,getIntervalsPercentileBoot);
	interFonseca = zero(0.005,100,machep,t,getIntervalsFonseca);
	interPlugIn = qweibull(alpha, betaMLE, etaMLE, 1, 0);
	REAL(MLEs_Data)[4+n+2*B+4] = interGPQ;
	REAL(MLEs_Data)[4+n+2*B+5] = interPerBoot;
	REAL(MLEs_Data)[4+n+2*B+6] = interFonseca;
	REAL(MLEs_Data)[4+n+2*B+7] = interPlugIn;
	
	UNPROTECT(1);
	return MLEs_Data;
}
