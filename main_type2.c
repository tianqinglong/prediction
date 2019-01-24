#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall_type2.h"
#include "brent.h"

int r=30, n=50;
double weiEta=1.5, weiBeta=1.5;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX], weightArray[ARRAY_MAX];

int main(){
	int i;
	double betaMLE, etaMLE, betaMLEBoot, etaMLEBoot;
	double t, machep, nt, dt;

	set_seed(time(NULL),500);
	simuDataType2(cenArray, realArray);

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.1,10,machep,t,deriEquation);

	etaMLE = 0;
	for(i=0;i<n;i++){
		etaMLE += pow(cenArray[i], betaMLE);
	}
	etaMLE = pow(etaMLE/r,1/betaMLE);

	printf("%f %f\n", betaMLE, etaMLE);

	getWeight(weightArray);
	/*
	for(i=1;i<n;i++){
		printf("%f\n", weightArray[i]);
	}*/
	betaMLEBoot = zero(0.1,10,machep,t,deriFRWB);
	nt=0;dt=0;
	for(i=0;i<n;i++){
		nt += weightArray[i]*pow(cenArray[i],betaMLEBoot);
		dt += ((i<r) ? weightArray[i] : 0);
	}
	etaMLEBoot = pow(nt/dt, 1/betaMLEBoot);
	printf("%f %f\n", betaMLEBoot, etaMLEBoot);

	return 0;
}