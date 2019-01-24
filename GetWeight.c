#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

extern int n;

void getWeight(double *weight){
	int i;
	double V=0,Y[n];

	for(i=0;i<n;i++){
		Y[i] = rexp(1);
		V += Y[i];
	}

	for(i=0;i<n;i++){
		weight[i] = Y[i]/V;
	}
}