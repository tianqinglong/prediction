#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "overall.h"

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

double deriEquation(double weibullBeta){
	int i;
	double equation, weibullEta=0;

	for(i=0;i<n;i++){
		weibullEta += pow(cenArray[i], weibullBeta);
	}
	weibullEta = pow(weibullEta/r,1/weibullBeta);

	equation = r/weibullBeta - r*log(weibullEta);
	for(i=0;i<r;i++){
		equation += log(cenArray[i]);
	}
	for(i=0;i<n;i++){
		equation -= pow(cenArray[i]/weibullEta, weibullBeta)*log(cenArray[i]/weibullEta);
	}

	return equation;
}

double deriEquationParaBoot(double shape){
	int i;
	double equation, scale=0;

	for(i=0;i<n;i++){
		scale += pow(cenBootSample[i], shape);
	}
	scale = pow(scale/r,1/shape);

	equation = r/shape - r*log(scale);
	for(i=0;i<r;i++){
		equation += log(cenBootSample[i]);
	}
	for(i=0;i<n;i++){
		equation -= pow(cenBootSample[i]/scale, shape)*log(cenBootSample[i]/scale);
	}

	return equation;
}

double deriFRWB(double weibullBeta){
	int i;
	double weibullEta, nt, dt, equation;

	nt=0;dt=0;
	for(i=0;i<n;i++){
		nt += weightArray[i]*pow(cenArray[i],weibullBeta);
		dt += ((i<r) ? weightArray[i] : 0);
	}
	weibullEta = pow(nt/dt, 1/weibullBeta);

	equation = dt/weibullBeta;
	for(i=0;i<r;i++){
		equation += weightArray[i]*log(cenArray[i]/weibullEta);
	}
	for(i=0;i<n;i++){
		equation -= weightArray[i]*pow(cenArray[i]/weibullEta, weibullBeta)*log(cenArray[i]/weibullEta);
	}

	return equation;
}
