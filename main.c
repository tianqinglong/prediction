#include <stdio.h>
#include <time.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "prediction.h"

double censor; // t_{c}

int main(){
	int i;
	int type, FRWB;
	int discrete;
	double Er, Pt, shape, scale, lower, upper, times;
	double *cptmp;
	double cp[N][3], cp1=0, cp2=0, cp3=0;

	type = 1; // Censoring Type: 1 or 2
	Er = 5;	// (Expected) number of failures
	Pt = 0.1; // (Expected) Fraction of failing
	shape = 1.5;
	scale = 1;
	FRWB = 0; // Use FRWB: 1 yes, 0 no
	lower = 0.05;
	upper = 0.95;
	times = 1.5; // In discrete case, t_{w} = times * t_{c}
	discrete = 1; // 1: discrete prediction, 0: continuous prediction

	set_seed(time(NULL),1818);

	if(discrete == 0)
	{
		for(i=0;i<N;i++)
		{
			printf("This is iteration %d.\n", i+1);
			cptmp = single_continous_iteration(type, Er, Pt, shape, scale, FRWB, lower, upper);
			cp[i][0] = cptmp[0];
			cp[i][1] = cptmp[1];
			cp[i][2] = cptmp[2];
		}

		printf("\t\tPlug-in \tGPQ \t\tPercentile\n");
		for(i=0;i<N;i++)
		{
			printf("Iteration %d:\t%f \t%f \t%f\n", i+1, cp[i][0], cp[i][1], cp[i][2]);

			cp1 += cp[i][0]; cp2 += cp[i][1]; cp3 += cp[i][2];
		}

		cp1 = cp1/N; cp2 = cp2/N; cp3 = cp3/N;

		printf("%f %f %f\n", cp1, cp2, cp3);

		printf("\a");
	}
	else if(discrete == 1)
	{
		for(i=0;i<N;i++){
			printf("This is iteration %d.\n", i+1);
			single_binom_iteration(type, Er, Pt, shape, scale, FRWB, lower, upper, times);
			printf("\n");
		}
		// printf("%f\n", censor);
	}

	return 0;
}
