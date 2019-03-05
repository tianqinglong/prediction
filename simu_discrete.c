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
	FILE *fp2 = NULL;
	FILE *fp1_2 = NULL;
	FILE *fp2_2 = NULL;

	fp1 = fopen("discete_simulation_type1.txt", "a+");
	fputs("PB\tPI\tGPQ\tFon\n", fp1);
	fp2 = fopen("index_discete_type1.txt", "a+");
	fputs("Type\tEr\tPt\tPt2\tBeta\tFRWB\tPB\tPI\tGPQ\tFon\n", fp2);

	fp1_2 = fopen("discete_simulation_type2.txt", "a+");
	fputs("PB\tPI\tGPQ\tFon\n", fp1_2);
	fp2_2 = fopen("index_discete_type2.txt", "a+");
	fputs("Type\tEr\tPt\tPt2\tBeta\tFRWB\tPB\tPI\tGPQ\tFon\n", fp2_2);

    int i;

    int type = 1, isFRWB = 0;
    double nextCenorP;
    double add = 0.1;
    double Er_array[2] = {5, 10}, Pt_array[2] = {0.05, 0.1};
    double lower = 0.025, upper = 0.975;
    double beta = 1, eta=1;
    double interval;
    double Er, Pt;

    int r_argc = 2;
	char *r_argv[] = { "R", "--silent" };
	Rf_initEmbeddedR(r_argc, r_argv);

	set_seed(1234,12);
	double *cptmp, *cptmp2;

	double cp1=0, cp2=0, cp3=0, cp4=0;
	double cp1_2 = 0, cp2_2=0, cp3_2=0, cp4_2=0;

	int ee, pp;
	for(ee = 0;ee < 2; ee++){

		Er = Er_array[ee];

		for(pp =0;pp <2;pp++){

			Pt = Pt_array[pp];

			// for(aa = 0 ; aa < 2 ; aa++){

				// add = add_array[aa];
				nextCenorP = Pt + add;

				// for(bb=0; bb < 2 ;bb++){

					// beta = beta_array[bb];

					cp1=0; cp2=0; cp3=0; cp4=0;


					interval = qweibull(nextCenorP, beta, eta, 1, 0) - qweibull(Pt, beta, eta, 1, 0);
					for(i=0;i<N;i++)
					{
						printf( "This is iteration %d\n", i+1 );
						cptmp = single_binom_iteration(type, Er, Pt, beta, eta, isFRWB, lower, upper, nextCenorP);
						cptmp2 = single_type2_binom( Er, Pt, beta, eta, isFRWB, lower, upper, interval);
						// temp = single_binom_iteration_calibration_type1( Er,  Pt,  beta,  eta,  lower,  upper,  nextCenorP);
						fprintf(fp1, "%f\t%f\t%f\t%f\n", cptmp[0], cptmp[1], cptmp[2], cptmp[3]);
						fprintf(fp1_2, "%f\t%f\t%f\t%f\n", cptmp2[0], cptmp2[1], cptmp2[2], cptmp2[3]);

						cp1 += cptmp[0];cp2 += cptmp[1];cp3 += cptmp[2];cp4 += cptmp[3];
						cp1_2 +=cptmp2[0]; cp2_2 += cptmp2[1];
						cp3_2 +=cptmp2[2]; cp4_2 += cptmp2[3];
					}

					cp1 = cp1/N; cp2 = cp2/N; cp3 = cp3/N; cp4 = cp4/N;

					cp1_2 = cp1_2/N; cp2_2 = cp2_2/N; cp3_2 = cp3_2/N; cp4_2 = cp4_2/N;
					fprintf(fp2, "%d\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n", type, (int) Er, Pt, nextCenorP, beta, isFRWB, cp1, cp2, cp3, cp4);
					fprintf(fp2_2, "%d\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n", type, (int) Er, Pt, nextCenorP, beta, isFRWB, cp1_2, cp2_2, cp3_2, cp4_2);

				// }
			// }
		}
	}
	

    Rf_endEmbeddedR(0);

    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );
    fclose(fp1);
    fclose(fp2);
    return 0;
}