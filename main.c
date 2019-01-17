#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#define n 20
#define r 5

void simuTypeTwo(double *numVect, int *cenVect);
double likTypeTwo(double *numVect, double paraEta, double paraBeta);
double typeTwoMaximizer(double paraBeta);

double weibullEta=3;
double weibullBeta=2;

double weibullArray[n];
int weibullCensor[n];

int main(){
	int i;

	set_seed(time(NULL),time(NULL)+1);
	simuTypeTwo(weibullArray, weibullCensor);

	for(i=1;i<=n;i++){
		printf("%f %d\n", weibullArray[i-1], weibullCensor[i-1]);
	}

	printf("%f %f\n", likTypeTwo(weibullArray, 3, 2), typeTwoMaximizer(2));

	return 0;
}

void simuTypeTwo(double *numVect, int *cenVect){
	int i;
	double u0, u, ord;

	u0=0;
	u = unif_rand();

	for (i=1; i<=n; i++){
		ord = 1 - (1-u0)*pow(1-u, 1.0/(n-(i-1)));
		u0 = ord;
		u = unif_rand();
		numVect[i-1] = weibullEta*pow(-log(1-ord),1/weibullBeta);
		cenVect[i-1] = !(i<=r);
	}
}

double likTypeTwo(double *numVect, double paraEta, double paraBeta){
	int i;
	double loglik, survProb;

	loglik = 0;

	for(i=1; i<=r; i++){
		loglik += dweibull(numVect[i-1], paraBeta, paraEta, TRUE);
	}
	survProb = pweibull(numVect[r-1], paraBeta, paraEta, FALSE, TRUE);
	loglik = loglik + (n-r)*survProb;

	return loglik;
}

double typeTwoMaximizer(double paraBeta){
	int i;
	double paraEta, sumOfArray, loglik, survProb;

	sumOfArray = 0;
	for(i=1;i<=r;i++){
		sumOfArray += pow(weibullArray[i-1],paraBeta);
	}
	sumOfArray += (n-r)*pow(weibullArray[r-1],paraBeta);
	paraEta = pow(sumOfArray/r,1/paraBeta);

	loglik = 0;

	for(i=1; i<=r; i++){
		loglik += dweibull(weibullArray[i-1], paraBeta, paraEta, TRUE);
	}
	survProb = pweibull(weibullArray[r-1], paraBeta, paraEta, FALSE, TRUE);
	loglik = loglik + (n-r)*survProb;

	return loglik;
}


