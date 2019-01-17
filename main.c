#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

void simuTypeTwo(double *numVect, int *cenVect);

int n=10;
int r=4;
double weibullEta=1;
double weibullBeta=1.5;

int main(){
	int i;
	double weibull[n];
	int censor[n];

	set_seed(time(NULL), 123);
	simuTypeTwo(weibull, censor);

	for(i=1;i<=n;i++){
		printf("%f %d\n", weibull[i-1], censor[i-1]);
	}
	
	return 0;
}

void simuTypeTwo(double *numVect, int *cenVect){
	int i;
	double u;

	for (i=1; i<=n; i++){
		u = unif_rand();
		numVect[i-1] = weibullEta*pow(-log(1-u),1/weibullBeta);
		cenVect[i-1] = i<=r;
	}
}