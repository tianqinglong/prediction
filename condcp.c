#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "prediction.h"

extern double condiprob(double shape, double scale, double times);
double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times);
int qualify_discrete(double betab, double etab, double beta_mle, double eta_mle, double times);
int qualify_continuous(double betab, double etab, double eta_mle, double beta_mle);

double *single_continous_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper)
// Compute the conditional coverage probability for one simulated dataset
{
	int i;

	int r,n;
	int r_db, n_db;
	double *data, *weight0, *MLEs, *mb, *wb, *db;
	double mbeta, meta; //local copies to record mles
	double betaB[B], etaB[B];

	static double cp[3];

	if(type == 1)
	{
		censor = eta*pow(-log(1-Pt), 1/beta);
		data = simulator(type, Er, Pt, beta, eta);
		r = data[0];
		n = data[1];
	}
	else
	{
		data = simulator(type, Er, Pt, beta, eta);
		r = data[0];
		n = data[1];
		censor = data[r+1];
	}

	printf("The number of failure is: %d\nThe sample size is: %d\n",r,n);
	
	// for(i=0;i<n;i++){
	// 	printf("%f\n", data[i+2]);
	// }

	weight0 = generateWeights(0, r, n); // Get the MLEs of the initial dataset, no need to do FRWB

	// for(i=0;i<r+1;i++){
	// 	printf("%f\n", weight0[i]);
	// }

	MLEs = findmle(data, weight0);

	// make a local copy
	mbeta = MLEs[0];
	meta = MLEs[1];
	printf("The MLEs are: (%f, %f)\n", mbeta, meta);

	if(FRWB == 1)
	{
		for(i=0;i<B;i++)
		{
			wb = generateWeights(1, r, n);
			// if(r==3 && n==60)
			// for(j=0;j<(r+1);j++){
			// 	printf("%f ", wb[j]);
			// }

			mb = findmle(data, wb);
			betaB[i] = mb[0];
			etaB[i] = mb[1];

			if(qualify_continuous(betaB[i], etaB[i], mbeta, meta)){
				printf("Abnormal Bootstrap Draws %f %f\n", betaB[i], etaB[i]);
				for (int j = 0; j < (r+1); ++j)
				{
					printf("%f %f\n", wb[j], data[j+2]);
				}
				getchar();
				i--;
			}
		}
	}
	else if(FRWB == 0)
	{
		if(type == 2)
		{
			for(i=0;i<B;i++)
			{
				db = simulator(type, Er, Pt, mbeta, meta);

				r_db = db[0];
				n_db = db[1];

				wb = generateWeights(0, r_db, n_db);

				mb = findmle(db, wb);
				betaB[i] = mb[0];
				etaB[i] = mb[1];

				if(qualify_continuous(betaB[i], etaB[i], mbeta, meta)){
					printf("Abnormal Bootstrap Draws %f %f\n", betaB[i], etaB[i]);
					i--;
				}
			}
		}
		else if(type == 1)
		{
			for(i=0;i<B;i++)
			{
				db = bootSimulator(Er, Pt, mbeta, meta);

				// for(i=0;i<r+1;i++){
				// 	printf("%f\n", weight1[i]);
				// }

				r_db = db[0];
				n_db = db[1];

				wb = generateWeights(0, r_db, n_db);

				mb = findmle(db, wb);
				betaB[i] = mb[0];
				etaB[i] = mb[1];

				if(qualify_continuous(betaB[i], etaB[i], mbeta, meta)){
					printf("Abnormal Bootstrap Draws %f %f\n", betaB[i], etaB[i]);
					i--;
					// int k;
					// for(k=0; k<n+2;k++)
					// {
					// 	printf("%f\n", db[k]);
					// }
				}
			}
		}
	}

	// printf("%f %f\n", MLEs[0], MLEs[1]);
	// for(i=0;i<B;i++)
	// {
	// 	printf("%f %f\n", betaB[i], etaB[i]);
	// }

	double ppil, ppiu; // plug-in
	double pcp, bcp, gcp; // condition coverage probablity for plug-in, percentile and gpq methods.

	// prediction interval for plug-in method
	ppil = qweibull(lower, mbeta, meta, 1, 0);
	ppiu = qweibull(upper, mbeta, meta, 1, 0);
	pcp = pweibull(ppiu, beta, eta, 1, 0) - pweibull(ppil, beta, eta, 1, 0);

	// printf("The plug-in prediction interval is: (%f,%f)\n", ppil, ppiu);
	// printf("The conditional probability of plug-in method is: %f\n\n", pcp);

	// prediction interval for GPQ method

	double *gpqinter; // gpq
	gpqinter = gpqinterval(betaB, etaB, lower, upper, mbeta, meta);
	gcp = pweibull(gpqinter[1], beta, eta, 1, 0) - pweibull(gpqinter[0], beta, eta, 1, 0);

	// printf("The gpq prediction interval is: (%f,%f)\n", gpil, gpiu);
	// printf("The conditional probability of gpq method is: %f\n\n", gcp);

	// prediction interval for Percentile Bootstrap

	double *bpinter; //bp
	bpinter = pbinterval(betaB, etaB, lower, upper);
	bcp = pweibull(bpinter[1], beta, eta, 1, 0) - pweibull(bpinter[0], beta, eta, 1, 0);

	// printf("The pb prediction interval is: (%f,%f)\n", bpil, bpiu);
	// printf("The conditional probability of pb method is: %f\n\n", bcp);

	cp[0] = pcp;
	cp[1] = gcp;
	cp[2] = bcp;

	//debugger(betaB, etaB);
	return cp;
}

double *single_binom_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper, double nextCen)
{
	int i;

	static double cp_binom[5];

	int r,n;
	double times;
	int ubinom_pb, lbinom_pb, ubinom_gpq, lbinom_gpq, ubinom_pi, lbinom_pi, ubinom_cal, lbinom_cal, ubinom_fon, lbinom_fon;
	double *data, *weight0, *MLEs, *mb, *wb, *db;
	double mbeta, meta; //local copies to record mles
	double betaB[B], etaB[B];

	double cp_pb, cp_gpq, cp_pi, cp_cali, cp_fon;

	//static double cp[3];

	if(type == 1)
	{
		censor = eta*pow(-log(1-Pt), 1/beta);
		data = simulator(type, Er, Pt, beta, eta);
		r = data[0];
		n = data[1];
	}
	else
	{
		data = simulator(type, Er, Pt, beta, eta);
		r = data[0];
		n = data[1];
		censor = data[r+1];
	}

	times = qweibull(nextCen, beta, eta, 1, 0)/censor;
	weight0 = generateWeights(0, r, n); // Get the MLEs of the initial dataset, no need to do FRWB

	// make a local copy
	double data1[n+2], weight1[n];
	for(i=0; i<(n+2); i++)
	{
		data1[i] = data[i];
		// printf("%f\n", data1[i]);
	}
	for(i=0; i<(r+1); i++)
	{
		weight1[i] = weight0[i];
	}

	MLEs = findmle(data1, weight1);

	// make a local copy
	mbeta = MLEs[0];
	meta = MLEs[1];
	printf("The MLEs are: (%f, %f)\n", mbeta, meta);

	if(FRWB == 1)
	{
		for(i=0;i<B;i++)
		{
			wb = generateWeights(1, r, n);
			mb = findmle(data1, wb);
			betaB[i] = mb[0];
			etaB[i] = mb[1];
			if(qualify_discrete(betaB[i], etaB[i], mbeta, meta, times))
			{
				i--;
			}
		}
	}
	else if(FRWB == 0)
	{
		if(type == 2)
		{
			for(i=0;i<B;i++)
			{
				db = simulator(type, Er, Pt, mbeta, meta);
				mb = findmle(db, weight1);
				betaB[i] = mb[0];
				etaB[i] = mb[1];

				if(qualify_discrete(betaB[i], etaB[i], mbeta, meta, times))
				{
					printf("Bootstrap Draws Abnormal");
					i--;
				}
			}
		}
		else if(type == 1)
		{
			for(i=0;i<B;i++)
			{
				db = bootSimulator(Er, Pt, mbeta, meta);
				mb = findmle(db, weight1);
				betaB[i] = mb[0];
				etaB[i] = mb[1];

				if(qualify_discrete(betaB[i], etaB[i], mbeta, meta, times))
				{
					printf("Bootstrap Draws Abnormal");
					i--;
				}
			}
		}
	}

	lbinom_pb = pbbinominterval(betaB, etaB, lower, times, n, r);
	ubinom_pb = pbbinominterval(betaB, etaB, upper, times, n, r);

	cp_pb = find_binom_prob(lbinom_pb, ubinom_pb, n-r, beta, eta, times);
	cp_binom[0] = cp_pb;

	printf("PB:[%d,%d] CP:%f\n", lbinom_pb, ubinom_pb, cp_pb);

	lbinom_gpq = gpqbinominterval(betaB, etaB, lower, times, mbeta, meta, n, r);
	ubinom_gpq = gpqbinominterval(betaB, etaB, upper, times, mbeta, meta, n, r);

	cp_gpq = find_binom_prob(lbinom_gpq, ubinom_gpq, n-r, beta, eta, times);
	cp_binom[1] = cp_gpq;

	printf("GPQ:[%d,%d] CP:%f\n", lbinom_gpq, ubinom_gpq, cp_gpq);

	double p_pi,tmp;

	p_pi = condiprob(mbeta, meta, times);
	tmp = qbinom(lower, n-r, p_pi,1,0);
	lbinom_pi = (tmp == 0) ? 0 : (tmp - 1);
	ubinom_pi = qbinom(upper, n-r, p_pi,1,0);

	cp_pi = find_binom_prob(lbinom_pi,ubinom_pi, n-r, beta, eta, times);
	cp_binom[2] = cp_pi;

	printf("PI:[%d,%d] CP:%f\n", lbinom_pi, ubinom_pi, cp_pi);

	double ucali, lcali;
	ucali = calibinominterval(betaB, etaB, upper, times, mbeta, meta, n, r);
	lcali = calibinominterval(betaB, etaB, lower, times, mbeta, meta, n, r);

	ubinom_cal = qbinom(ucali, n-r, p_pi, 1, 0);
	lbinom_cal = qbinom(lcali, n-r, p_pi, 1, 0);
	if(lbinom_cal > 0){
		lbinom_cal--;
	}
	cp_cali = find_binom_prob(lbinom_cal,ubinom_cal, n-r, beta, eta, times);
	cp_binom[3] = cp_cali;

	printf("CALI:[%d,%d] CP:%f\n", lbinom_cal, ubinom_cal, cp_cali);

	ubinom_fon = fonsecabinominterval(betaB, etaB, upper, times, mbeta, meta, n, r);
	lbinom_fon = fonsecabinominterval(betaB, etaB, lower, times, mbeta, meta, n, r);

	cp_fon = find_binom_prob(lbinom_fon, ubinom_fon, n-r, beta, eta, times);
	cp_binom[4] = cp_fon;

	printf("FON:[%d,%d] CP:%f\n\n", lbinom_fon, ubinom_fon, cp_fon);

	return cp_binom;
}

double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times)
{
	double value, prob;

	prob = condiprob(beta, eta, times);
	value = pbinom(upper, n, prob, 1, 0) - pbinom(lower, n, prob, 1, 0) + dbinom(lower, n, prob, 0);

	return value;
}

int qualify_discrete(double betab, double etab, double beta_mle, double eta_mle, double times)
{
	int flag1, flag2;

	double binomp;
	binomp = condiprob(betab, etab, times);
	flag1 = isnan(binomp);

	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	mub = log(etab);
	sigmab = 1/betab;

	betagpq = 1/(sigmah/sigmab*sigmah);
	etagpq = exp(muh + (muh-mub)/sigmab*sigmah);
	
	binomp = condiprob(betagpq, etagpq, times);

	flag2 = isnan(binomp);

	return (flag1||flag2);
}

int qualify_continuous(double betab, double etab, double eta_mle, double beta_mle)
{
	int flag2, flag3;
	double cdf0, cdf1;
	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	cdf0 = pweibull(0.01, betab, etab, 1, 0);
	flag2 = isnan(cdf0);

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	mub = log(etab);
	sigmab = 1/betab;

	betagpq = 1/(sigmah/sigmab*sigmah);
	etagpq = exp(muh + (muh-mub)/sigmab*sigmah);
	
	cdf1 = pweibull(0.01, betagpq, etagpq, 1 ,0);

	flag3 = isnan(cdf1);

	return (flag3||flag2);
}
