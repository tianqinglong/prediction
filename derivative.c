#include <stdio.h>
#include <math.h>

extern int r, n;
extern double cenArray[];

double deriEquation1(double weibullBeta){
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