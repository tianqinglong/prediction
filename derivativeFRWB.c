#include <stdio.h>
#include <math.h>

extern int r,n;
extern double weightArray[], cenArray[];

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