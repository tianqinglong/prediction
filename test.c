#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "type2censor.h"
#include "brent.h"

int n=50, r=30;
double weiEta=1.5, weiBeta=1.5;
double cenArray[ARRAY_MAX], realArray[ARRAY_MAX];

int main(){
	int i;
	double betaMLE, etaMLE;
	double t, machep;

	set_seed(time(NULL), 580580);
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

	return 0;
}