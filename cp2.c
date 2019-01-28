#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <R.h>
#include <Rinternals.h>

#include "overall.h"
#include "brent.h"

int r, n;
double Er, Pt;
double weiEta, weiBeta, alpha;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX], weightArray[ARRAY_MAX];
double cenBootSample[ARRAY_MAX], comBootSample[ARRAY_MAX];
double etaMLE, betaMLE;
double betaBoot[B], etaBoot[B];

SEXP transferDataToR(SEXP Rseed, SEXP lower, SEXP upper, SEXP Rr, SEXP Rn, SEXP Reta, SEXP Rbeta, SEXP RBoot){
	int i,j;
	int boottype;
    
    Er = r;
    Pt = (float)r/n;
    
	int Seed = asInteger(Rseed);
	r = asInteger(Rr), n = asInteger(Rn);
	double lAlpha = asReal(lower), uAlpha = asReal(upper);
	weiEta = asReal(Reta), weiBeta = asReal(Rbeta);
	boottype = asInteger(RBoot);

	double betaMLEBoot, etaMLEBoot;
	double interGPQ, interPerBoot, interFonseca, interPlugIn;
	double t, machep, nt, dt;

	SEXP MLEs_Data = PROTECT(allocVector(REALSXP,(4+n+2*B+8)));

	set_seed(time(NULL)/Seed,Seed*time(NULL));
	simuDataType2(cenArray, realArray, weiBeta, weiEta);

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.005,10,machep,t,deriEquation);

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
	
	if(boottype){
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
		for(j=0;j<B;j++){
			simuDataType2(cenBootSample, comBootSample, betaMLE, etaMLE);

			betaMLEBoot = zero(0.01,100,machep,t,deriEquationParaBoot);
		
			etaMLEBoot = 0;
			for(i=0;i<n;i++){
				etaMLEBoot += pow(cenBootSample[i], betaMLEBoot);
			}
			etaMLEBoot = pow(etaMLEBoot/r,1/betaMLEBoot);

			betaBoot[j] = betaMLEBoot;
			etaBoot[j] = etaMLEBoot;

			REAL(MLEs_Data)[4+n+j] = betaMLEBoot;
			REAL(MLEs_Data)[4+n+B+j] = etaMLEBoot;
		}
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
