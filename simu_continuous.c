#include <stdio.h>
#include <time.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <Rembedded.h>

#include "prediction.h"

double censor; // t_{c}

int main()
{
	FILE *fp1 = NULL;
	fp1 = fopen("continuous_simulation.txt", "a+");

	fputs("Type\tEr\tPt\tBeta\tPlug\tGPQ\tPB\n", fp1);

	printf("Type\tEr\tPt\tBeta\tPlug\tGPQ\tPB\n");

	int i,j,k,m;
	int type = 1;

	double *cptmp;
	double cp1=0, cp2=0, cp3=0;

	double eta=1;
	double beta, Er, Pf;
	double beta_array[1]={1}, Er_array[3]={5, 10, 15}, Pf_array[3]={0.05 ,0.1 ,0.2};

	set_seed(12, 1818);

	int r_argc = 2;
	char *r_argv[] = { "R", "--silent" };
    Rf_initEmbeddedR(r_argc, r_argv);

    for(i=0;i<3;i++)
    {
    	beta = beta_array[i];

    	for(j=0;j<3;j++)
    	{
    		Er = Er_array[j];

    		for(k=0;k<3;k++)
    		{
    			Pf = Pf_array[k];

    			for(m=0;m<N;m++)
    			{
    				cptmp = single_continous_iteration(type, Er, Pf, beta, eta, 0, 0.025, 0.975);
    				cp1 += cptmp[0];
    				cp2 += cptmp[1];
    				cp3 += cptmp[2];
    			}

    			cp1 = cp1/N;
    			cp2 = cp2/N;
    			cp3 = cp3/N;

    			printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n", type, (int) Er, Pf, beta, cp1, cp2, cp3);
    			fprintf(fp1, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", type, (int) Er, Pf, beta, cp1, cp2, cp3);

    			cp1 = 0; cp2 = 0; cp3 = 0;
    		}
    	}
    }

    Rf_endEmbeddedR(0);

    fclose(fp1);
	return 0;
}