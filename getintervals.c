#include <stdio.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall.h"
#include "brent.h"

double getIntervalsGPQ(double quantile){
	int i;
	double muBoot, sigmaBoot, muHat, sigmaHat, muGPQ, sigmaGPQ;
	double equation=0, betaGPQ, etaGPQ;

	for(i=0;i<B;i++){
		muBoot = log(etaBoot[i]);
		sigmaBoot = 1/betaBoot[i];

		muHat = log(etaMLE);
		sigmaHat = 1/betaMLE;

		muGPQ = muHat + (muHat - muBoot)/sigmaBoot*sigmaHat;
		sigmaGPQ = sigmaHat/sigmaBoot*sigmaHat;

		betaGPQ = 1/sigmaGPQ; etaGPQ = exp(muGPQ);

		equation += pweibull(quantile, betaGPQ, etaGPQ, 1, 0);
	}
	equation = equation/B - alpha;

	return equation;
}

double getIntervalsFonseca(double quantile){
	int i;
	double equation=0;

	for(i=0;i<B;i++){
		equation += pweibull(qweibull(pweibull(quantile, betaMLE, etaMLE, 1, 0)
			,betaBoot[i],etaBoot[i],1,0),betaMLE, etaMLE, 1, 0);
	}
	equation = equation/B - alpha;

	return equation;
}

double getIntervalsPercentileBoot(double quantile){
	int i;
	double equation=0;

	for(i=0;i<B;i++){
		equation += pweibull(quantile, betaBoot[i], etaBoot[i], 1, 0);
	}
	equation = equation/B - alpha;

	return equation;
}
