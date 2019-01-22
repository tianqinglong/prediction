#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "brent.h"
#include "overall_type1.h"

int r,n;
double Pt=0.25, Er=10;
double weiEta=1.5, weiBeta=1.5;

double cenArray[ARRAY_MAX], comArray[ARRAY_MAX];

int main(){
	int i;
	double t, machep;
	double betaMLE, etaMLE;

	set_seed(time(NULL),500);

	r=0;
	while(r < 2){
		simuDataType1(cenArray, comArray);
	}

	for(i=0;i<n;i++){
		printf("%f %f\n", cenArray[i], comArray[i]);
	}
	printf("%d %d\n", r, n);

	machep = r8_epsilon();
	t = machep;

	betaMLE = zero(0.1,10,machep,t,deriEquation1);

	etaMLE = 0;
	for(i=0;i<n;i++){
		etaMLE += pow(cenArray[i], betaMLE);
	}
	etaMLE = pow(etaMLE/r,1/betaMLE);

	printf("%f %f\n", betaMLE, etaMLE);
	return 0;
}