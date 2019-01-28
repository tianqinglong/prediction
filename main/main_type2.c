#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <unistd.h>

#include "overall_type2.h"
#include "brent.h"

int r=30, n=100;
double weiEta=1.5, weiBeta=1.5, alpha = 0.05;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX], weightArray[ARRAY_MAX];
double etaMLE, betaMLE;
double betaBoot[B], etaBoot[B];

//gcc main_type2.c brent.c simulator_type2.c derivative.c derivativeFRWB.c getintervals.c GetWeight.c -Wall -lRmath -lm -o "main2.out"

int main(){
	int i,j;
	double betaMLEBoot, etaMLEBoot;
	double interGPQ, interPerBoot, interFonseca;
	double t, machep, nt, dt;

	set_seed(time(NULL),getpid());
	simuDataType2(cenArray, realArray, weiBeta, weiEta);
/*
	printf("%d %d\n", r, n);
*/
	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.1,10,machep,t,deriEquation);

	etaMLE = 0;
	for(i=0;i<n;i++){
		etaMLE += pow(cenArray[i], betaMLE);
	}
	etaMLE = pow(etaMLE/r,1/betaMLE);
/*
	printf("The MLEs are:\n%f %f\n", betaMLE, etaMLE);

	printf("The MLEs of Bootstrap Samples are:\n");
*/
	for(j=0;j<B;j++){
		getWeight(weightArray);
		betaMLEBoot = zero(0.005,100,machep,t,deriFRWB);
		nt=0;dt=0;
		for(i=0;i<n;i++){
			nt += weightArray[i]*pow(cenArray[i],betaMLEBoot);
			dt += ((i<r) ? weightArray[i] : 0);
		}
		etaMLEBoot = pow(nt/dt, 1/betaMLEBoot);
		betaBoot[j] = betaMLEBoot;
		etaBoot[j] = etaMLEBoot;
	}

	interGPQ = zero(0.01,100,machep,t,getIntervalsGPQ);
	interPerBoot = zero(0.01,100,machep,t,getIntervalsPercentileBoot);
	interFonseca = zero(0.01,100,machep,t,getIntervalsFonseca);
	printf("%f %f %f\n", interGPQ, interPerBoot, interFonseca);

	return 0;
}