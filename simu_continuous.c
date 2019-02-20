#include <stdio.h>
#include <time.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <Rembedded.h>

#include "prediction.h"

double censor; // t_{c}

int main()
{
    clock_t start, end;
    double duration;

    start = clock();
	FILE *fp1 = NULL;
	fp1 = fopen("continuous_simulation.txt", "a+");

    FILE *fp2 = NULL;
    fp2 = fopen("cp_continuous.txt", "a+");
    fprintf(fp2,"Plug GPQ PB\n");

	fputs("Type\tEr\tPt\tBeta\tPlug\tGPQ\tPB\n", fp1);

	printf("Type\tEr\tPt\tBeta\tPlug\tGPQ\tPB\n");

	int i,j,k,m;
	int type = 1;
    int isFRWB = 1;
    double lower = 0.025, upper = 0.975;

	double *cptmp;
	double cp1=0, cp2=0, cp3=0, cp_array[N][3];

	double eta=1;
	double beta, Er, Pf;

    int lenbeta = 1, lenEr = 4, lenPf = 3;
	// double beta_array[lenbeta]={0.8, 1, 1.5, 3}, Er_array[lenEr]={5, 10, 15}, Pf_array[lenPf]={0.05 ,0.1 ,0.2};

    double beta_array[1]={1}, Er_array[4]={3, 5, 10, 15}, Pf_array[3]={0.05 ,0.1, 0.2};

	set_seed(12, 1818);

	int r_argc = 2;
	char *r_argv[] = { "R", "--silent" };
    Rf_initEmbeddedR(r_argc, r_argv);

    for(i=0;i<lenbeta;i++)
    {
    	beta = beta_array[i];

    	for(j=0;j<lenEr;j++)
    	{
    		Er = Er_array[j];

    		for(k=0;k<lenPf;k++)
    		{
    			Pf = Pf_array[k];

    			for(m=0;m<N;m++)
    			{
    				cptmp = single_continous_iteration(type, Er, Pf, beta, eta, isFRWB, lower, upper);
    				cp1 += cptmp[0];
    				cp2 += cptmp[1];
    				cp3 += cptmp[2];
                    cp_array[m][0] = cptmp[0];
                    cp_array[m][1] = cptmp[1];
                    cp_array[m][2] = cptmp[2];
    			}

                for(m=0;m<N;m++)
                {
                    fprintf(fp2, "%f %f %f\n", cp_array[m][0], cp_array[m][1], cp_array[m][2]);
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

    fclose(fp2);
    fclose(fp1);

    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );
	return 0;
}