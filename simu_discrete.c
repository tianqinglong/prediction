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
	fputs("Type\tFRWB\tEr\tPt\tPB\tPI\tGPQ\tFon\tNomi\n", fp1);
	fp2 = fopen("index_discete_type1.txt", "a+");
	fputs("Type\tEr\tPt\tPt2\tBeta\tFRWB\tPB\tPI\tGPQ\tFon\tNomi\n", fp2);

	fp1_2 = fopen("discete_simulation_type2.txt", "a+");
	fputs("Type\tFRWB\tEr\tPt\tPB\tPI\tGPQ\tFon\tNomi\n", fp1_2);
	fp2_2 = fopen("index_discete_type2.txt", "a+");
	fputs("Type\tEr\tPt\tPt2\tBeta\tFRWB\tPB\tPI\tGPQ\tFon\tNomi\n", fp2_2);

    int i;

    int type = 1, isFRWB = 0;
    double nextCenorP;
    double add = 0.1;
    double Er_array[7] = {5, 10, 15, 20, 25, 35, 50, 60}, Pt_array[2] = {0.05, 0.1};
    double lower = 0, upper1 = 0.95, upper2 = 0.9, upper3 = 0.05, upper4 = 0.1;
    double beta = 1, eta=1;
    double interval;
    double Er, Pt;

    int r_argc = 2;
	char *r_argv[] = { "R", "--silent" };
	Rf_initEmbeddedR(r_argc, r_argv);

	set_seed(1234,12);
	double *cptmp, *cptmp2;

	double cp1_11 = 0, cp1_21 = 0, cp1_31 = 0, cp1_41 = 0, cp1_12 = 0, cp1_22 = 0, cp1_32 = 0, cp1_42 = 0, cp1_13 = 0, cp1_23 = 0, cp1_33 = 0, cp1_43 = 0, cp1_14 = 0, cp1_24 = 0, cp1_34 = 0, cp1_44 = 0;
	double cp2_11 = 0, cp2_21 = 0, cp2_31 = 0, cp2_41 = 0, cp2_12 = 0, cp2_22 = 0, cp2_32 = 0, cp2_42 = 0, cp2_13 = 0, cp2_23 = 0, cp2_33 = 0, cp2_43 = 0, cp2_14 = 0, cp2_24 = 0, cp2_34 = 0, cp2_44 = 0;

	int ee, pp;
	for(ee = 0 ; ee < 7 ; ee++){

		Er = Er_array[ee];

		for(pp =0;pp <2;pp++){

			Pt = Pt_array[pp];

			nextCenorP = Pt + add;

			cp1_11 = 0; cp1_21 = 0; cp1_31 = 0; cp1_41 = 0; cp1_12 = 0; cp1_22 = 0; cp1_32 = 0; cp1_42 = 0; cp1_13 = 0; cp1_23 = 0; cp1_33 = 0; cp1_43 = 0; cp1_14 = 0; cp1_24 = 0; cp1_34 = 0; cp1_44 = 0;
			cp2_11 = 0; cp2_21 = 0; cp2_31 = 0; cp2_41 = 0; cp2_12 = 0; cp2_22 = 0; cp2_32 = 0; cp2_42 = 0; cp2_13 = 0; cp2_23 = 0; cp2_33 = 0; cp2_43 = 0; cp2_14 = 0; cp2_24 = 0; cp2_34 = 0; cp2_44 = 0;


			interval = qweibull(nextCenorP, beta, eta, 1, 0) - qweibull(Pt, beta, eta, 1, 0);
			for(i=0;i<N;i++)
			{
				printf( "This is iteration %d\n", i+1 );
				cptmp = single_binom_iteration(type, Er, Pt, beta, eta, isFRWB, lower, upper1, upper2, upper3, upper4, nextCenorP);
				fprintf(fp2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", type, Er, Pt, nextCenorP, beta, isFRWB, cptmp[0], cptmp[1], cptmp[2], cptmp[3], upper1);
				fprintf(fp2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", type, Er, Pt, nextCenorP, beta, isFRWB, cptmp[4], cptmp[5], cptmp[6], cptmp[7], upper2);
				fprintf(fp2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", type, Er, Pt, nextCenorP, beta, isFRWB, cptmp[8], cptmp[9], cptmp[10], cptmp[11], upper3);
				fprintf(fp2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", type, Er, Pt, nextCenorP, beta, isFRWB, cptmp[12], cptmp[13], cptmp[14], cptmp[15], upper4);

				cp1_11 += cptmp[0]; cp1_21 += cptmp[1]; cp1_31 += cptmp[2]; cp1_41 += cptmp[3];
				cp1_12 += cptmp[4]; cp1_22 += cptmp[5]; cp1_32 += cptmp[6]; cp1_42 += cptmp[7];
				cp1_13 += cptmp[8]; cp1_23 += cptmp[9]; cp1_33 += cptmp[10]; cp1_43 += cptmp[11];
				cp1_14 += cptmp[12]; cp1_24 += cptmp[13]; cp1_34 += cptmp[14]; cp1_44 += cptmp[15];


				cptmp2 = single_type2_binom( Er, Pt, beta, eta, isFRWB, lower, upper1, upper2, upper3, upper4, interval);
				fprintf(fp2_2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", 2, Er, Pt, nextCenorP, beta, isFRWB, cptmp2[0], cptmp2[1], cptmp2[2], cptmp2[3], upper1);
				fprintf(fp2_2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", 2, Er, Pt, nextCenorP, beta, isFRWB, cptmp2[4], cptmp2[5], cptmp2[6], cptmp2[7], upper2);
				fprintf(fp2_2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", 2, Er, Pt, nextCenorP, beta, isFRWB, cptmp2[8], cptmp2[9], cptmp2[10], cptmp2[11], upper3);
				fprintf(fp2_2, "%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n", 2, Er, Pt, nextCenorP, beta, isFRWB, cptmp2[12], cptmp2[13], cptmp2[14], cptmp2[15], upper4);

				cp2_11 += cptmp2[0]; cp2_21 += cptmp2[1]; cp2_31 += cptmp2[2]; cp2_41 += cptmp2[3];
				cp2_12 += cptmp2[4]; cp2_22 += cptmp2[5]; cp2_32 += cptmp2[6]; cp2_42 += cptmp2[7];
				cp2_13 += cptmp2[8]; cp2_23 += cptmp2[9]; cp2_33 += cptmp2[10]; cp2_43 += cptmp2[11];
				cp2_14 += cptmp2[12]; cp2_24 += cptmp2[13]; cp2_34 += cptmp2[14]; cp2_44 += cptmp2[15];
			}

			fprintf(fp1, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1, isFRWB,Er, Pt, cp1_11/N, cp1_21/N, cp1_31/N, cp1_41/N, upper1);
			fprintf(fp1, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1, isFRWB,Er, Pt, cp1_12/N, cp1_22/N, cp1_32/N, cp1_42/N, upper2);
			fprintf(fp1, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1, isFRWB,Er, Pt, cp1_13/N, cp1_23/N, cp1_33/N, cp1_43/N, upper3);
			fprintf(fp1, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1, isFRWB,Er, Pt, cp1_14/N, cp1_24/N, cp1_34/N, cp1_44/N, upper4);

			fprintf(fp1_2, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 2, isFRWB, Er, Pt, cp2_11/N, cp2_21/N, cp2_31/N, cp2_41/N, upper1);
			fprintf(fp1_2, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 2, isFRWB, Er, Pt, cp2_12/N, cp2_22/N, cp2_32/N, cp2_42/N, upper2);
			fprintf(fp1_2, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 2, isFRWB, Er, Pt, cp2_13/N, cp2_23/N, cp2_33/N, cp2_43/N, upper3);
			fprintf(fp1_2, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 2, isFRWB, Er, Pt, cp2_14/N, cp2_24/N, cp2_34/N, cp2_44/N, upper4);
		}

	}

    Rf_endEmbeddedR(0);

    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );
    fclose(fp1);
    fclose(fp2);
    fclose(fp1_2);
    fclose(fp2_2);
    return 0;
}