#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "prediction.h"

extern double condiprob_type2( double shape, double scale, double censoring, double duration);
extern double condiprob(double shape, double scale, double times);

double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times);
double find_binom_prob_type2( double lower, double upper, int n, double beta, double eta, double start, double duration);

int qualify_discrete(double betab, double etab, double beta_mle, double eta_mle, double times);
int qualify_continuous(double betab, double etab, double eta_mle, double beta_mle);
int qualify_discrete_type2(double betab, double etab, double beta_mle, double eta_mle, double start, double duration);

double get_real_plug_in_cover(int r, int n, double m_beta, double m_eta, double b_beta[], double b_eta[], double times, double lower, double upper);

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

	weight0 = generateWeights(0, r, n); // Get the MLEs of the initial dataset, no need to do FRWB
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
				}
			}
		}
	}

	// plug-in
	double ppil, ppiu;
	double pcp, bcp, gcp; // condition coverage probablity for plug-in, percentile and gpq methods.

	// prediction interval for plug-in method
	ppil = qweibull(lower, mbeta, meta, 1, 0);
	ppiu = qweibull(upper, mbeta, meta, 1, 0);
	pcp = pweibull(ppiu, beta, eta, 1, 0) - pweibull(ppil, beta, eta, 1, 0);

	// prediction interval for GPQ method

	double *gpqinter; // gpq
	gpqinter = gpqinterval(betaB, etaB, lower, upper, mbeta, meta);
	gcp = pweibull(gpqinter[1], beta, eta, 1, 0) - pweibull(gpqinter[0], beta, eta, 1, 0);

	// prediction interval for Percentile Bootstrap

	double *bpinter; //bp
	bpinter = pbinterval(betaB, etaB, lower, upper);
	bcp = pweibull(bpinter[1], beta, eta, 1, 0) - pweibull(bpinter[0], beta, eta, 1, 0);

	cp[0] = pcp;
	cp[1] = gcp;
	cp[2] = bcp;

	return cp;
}

double *single_binom_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper, double nextCen)
{
	int i;

	int r,n;
	int r_db, n_db;
	double *data, *weight0, *MLEs, *mb, *wb, *db;
	double mbeta, meta; //local copies to record mles
	double betaB[B], etaB[B];

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

	weight0 = generateWeights(0, r, n); // Get the MLEs of the initial dataset, no need to do FRWB
	MLEs = findmle(data, weight0);

	// make a local copy
	mbeta = MLEs[0];
	meta = MLEs[1];
	printf("The MLEs are: (%f, %f)\n", mbeta, meta);

	double times;
	times = qweibull(nextCen, beta, eta, 1, 0) / censor;

	if(FRWB == 1)
	{
		for(i=0;i<B;i++)
		{
			wb = generateWeights(1, r, n);

			mb = findmle(data, wb);
			betaB[i] = mb[0];
			etaB[i] = mb[1];

			if( qualify_discrete( mb[0], mb[1], mbeta, meta, times ) ){
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

			}
		}
		else if(type == 1)
		{
			for(i=0;i<B;i++)
			{
				db = bootSimulator(Er, Pt, mbeta, meta);

				r_db = db[0];
				n_db = db[1];

				wb = generateWeights(0, r_db, n_db);

				mb = findmle(db, wb);
				betaB[i] = mb[0];
				etaB[i] = mb[1];

				if( qualify_discrete( mb[0], mb[1], mbeta, meta, times ) ){
					i--;
				}

			}
		}
	}

	static double cp_binom[4];
	
	// printf("The time is %f\n", times);

// percentile bootstrap
	int *pb_binom_interval;
	pb_binom_interval = pbbinominterval(betaB, etaB, lower, upper, times, n, r);
	printf("PB: [%d,%d]\n", pb_binom_interval[0], pb_binom_interval[1]);
	cp_binom[0] = find_binom_prob(pb_binom_interval[0], pb_binom_interval[1], n-r, beta, eta, times);

// plug-in method
	int *plug_in_binom_interval;
	plug_in_binom_interval = plug_in_binominterval(mbeta, meta, lower, upper, times, n, r);
	printf("PI: [%d,%d]\n", plug_in_binom_interval[0], plug_in_binom_interval[1]);
	cp_binom[1] = find_binom_prob(plug_in_binom_interval[0],plug_in_binom_interval[1], n-r, beta, eta, times);

// GPQ method
	int *gpq_binom_interval;
	gpq_binom_interval = gpqbinominterval(betaB, etaB, lower, upper, times, mbeta, meta, n, r);
	printf("GPQ: [%d,%d]\n", gpq_binom_interval[0], gpq_binom_interval[1]);
	cp_binom[2] = find_binom_prob(gpq_binom_interval[0],gpq_binom_interval[1], n-r, beta, eta, times);

// Fonseca method
	int *fonseca_binom_interval;
	fonseca_binom_interval = fonsecabinominterval(betaB, etaB, lower, upper, times, mbeta, meta, n, r);
	printf("Fsc: [%d,%d]\n", fonseca_binom_interval[0], fonseca_binom_interval[1]);
	cp_binom[3] = find_binom_prob(gpq_binom_interval[0],gpq_binom_interval[1], n-r, beta, eta, times);

	return cp_binom;

}

// Find discete distribution for type 2 censoring
double *single_type2_binom(double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper, double interval)
{
	int type = 2;
	int i, r, n;

	double *data, *weight0, *MLEs;
	double beta_mle, eta_mle;

	data = simulator( type, Er, Pt, beta, eta );
	r = data[0];
	n = data[1];

	double data_local[n];
	// make a local copy
	for(i=0;i<(n+2);i++)
	{
		data_local[i] = data[i];
	}

	weight0 = generateWeights(0, r, n);

	MLEs = findmle(data, weight0);
	beta_mle = MLEs[0];
	eta_mle = MLEs[1];

	// Doing Bootstrap
	double *data_temp, *weight_temp, *mle_temp, beta_boot[B], eta_boot[B], censor_time[B];
	if(FRWB == 0)
	{
		weight_temp = generateWeights(0, r, n);

		for(i = 0; i < B; i++)
		{
			data_temp = simulator( type, Er, Pt, beta_mle, eta_mle );
			mle_temp = findmle( data_temp, weight_temp );

			beta_boot[i] = mle_temp[0];
			eta_boot[i] = mle_temp[1];
			censor_time[i] = data_temp[r+2];

			if( qualify_discrete_type2(mle_temp[0], mle_temp[1], beta_mle, eta_mle, censor_time[i], interval) )
			{
				i--;
			}
		}

	}
	else
	{
		for(i = 0; i < B; i++)
		{

			weight_temp = generateWeights(1, r, n);
			mle_temp = findmle(data_local, weight_temp);

			beta_boot[i] = mle_temp[0];
			eta_boot[i] = mle_temp[1];
			censor_time[i] = data_local[r+2];

			if( qualify_discrete_type2(mle_temp[0], mle_temp[1], beta_mle, eta_mle, censor_time[i], interval) )
			{
				i--;
			}
		}
	}

	static double cp_binom[4];

	// 2. Plug-in Method
	int *plug_in_binom_interval;
	printf("This is for Type 2:\n");
	plug_in_binom_interval = plug_in_binominterval_type2(beta_mle, eta_mle, lower, upper, data_local[r+2], interval, n, r);
	printf("PI: [%d,%d]\n", plug_in_binom_interval[0], plug_in_binom_interval[1]);
	cp_binom[1] = find_binom_prob_type2(plug_in_binom_interval[0], plug_in_binom_interval[1], n-r, beta, eta, data_local[r+2], interval);

	// 1. Percentile Bootstrap
	int *pb_binom_interval;
	pb_binom_interval = pbbinominterval_type2( beta_boot , eta_boot , censor_time, interval , lower, upper, n, r);
	printf("PB: [%d,%d]\n", pb_binom_interval[0], pb_binom_interval[1]);
	cp_binom[0] = find_binom_prob_type2(pb_binom_interval[0], pb_binom_interval[1], n-r, beta, eta, data_local[r+2], interval);

	// 3. GPQ
	int *gpq_binom_interval;
	gpq_binom_interval = gpqbinominterval_type2( beta_boot, eta_boot, censor_time, lower, upper, interval, beta_mle, eta_mle, n, r);
	printf("GPQ: [%d,%d]\n", gpq_binom_interval[0], gpq_binom_interval[1]);
	cp_binom[2] = find_binom_prob_type2(gpq_binom_interval[0], gpq_binom_interval[1], n-r, beta, eta, data_local[r+2], interval);

	// 4. Fonseca
	int *fonseca_binom_interval;
	fonseca_binom_interval = fonsecabinominterval_type2(beta_boot, eta_boot, censor_time, lower, upper, interval, beta_mle, eta_mle, n, r, data_local[r+2]);
	printf("Fon: [%d,%d]\n", fonseca_binom_interval[0], fonseca_binom_interval[1]);
	cp_binom[3] = find_binom_prob_type2(fonseca_binom_interval[0], fonseca_binom_interval[1], n-r, beta, eta, data_local[r+2], interval);

	return cp_binom;
}

double single_binom_iteration_calibration_type1(double Er, double Pt, double beta, double eta, double lower, double upper, double nextCen)
{
	int i;
	int r, n, br, bn;
	double *weight0, *MLEs, *bdata, *data;
	double coverage_probability=1;
	double lcl, ucl;
	double cp_binom;
	double mbeta, meta, bmbeta[B], bmeta[B];

	censor = qweibull(Pt, beta, eta, 1, 0);

	double times = qweibull(nextCen, beta, eta, 1, 0) / censor;

	data = simulator(1 , Er, Pt, beta, eta);
	r = data[0];
	n = data[1];

	int r_array[B];

	weight0 = generateWeights(0, r, n); // Get the MLEs of the initial dataset, no need to do FRWB
	MLEs = findmle(data, weight0);

	mbeta = MLEs[0];
	meta = MLEs[1];

	for(i = 0; i < B; i++)
	{
		bdata = bootSimulator( Er, Pt, mbeta, meta);
		br = bdata[0];
		bn = bdata[1];

		weight0 = generateWeights(0, br, bn);
		MLEs = findmle(bdata, weight0);

		bmbeta[i] = MLEs[0];
		bmeta[i] = MLEs[1];
		r_array[i] = br;
	}

	lcl = 0;
	ucl = 1;

	while(coverage_probability > (upper - lower))
	{
		lcl += 0.001;
		ucl -= 0.001;
		coverage_probability = get_real_plug_in_cover(r, n, mbeta, meta, bmbeta, bmeta, times, lcl, ucl);
		printf("%f\n", coverage_probability);
	}

	int *plug_in_binom_interval;
	plug_in_binom_interval = plug_in_binominterval(mbeta, meta, lcl, ucl, times, n, r);
	printf("Cali: [%d,%d]\n", plug_in_binom_interval[0], plug_in_binom_interval[1]);
	printf("Calibrated: %f %f %f\n", lcl, ucl, coverage_probability);
	cp_binom = find_binom_prob(plug_in_binom_interval[0],plug_in_binom_interval[1], n-r, beta, eta, times);

	return cp_binom;
}

double get_real_plug_in_cover(int r, int n, double m_beta, double m_eta, double b_beta[], double b_eta[], double times, double lower, double upper)
{
	int i;
	int *plug_in_interval;
	double cp = 0;
	for(i=0;i<B;i++)
	{
		plug_in_interval = plug_in_binominterval(b_beta[i], b_eta[i], lower, upper, times, n, r);
		// prob = condiprob(b_beta[i], b_eta[i], times);
		// printf("%f %f %f %f %f %d %d %f\n", b_beta[i], b_eta[i], lower, upper, times, n, r_array[i], prob);
		// printf("%d %d %d %f %f\n", plug_in_interval[0], plug_in_interval[1], n-r, lower, upper);
		cp += find_binom_prob(plug_in_interval[0],plug_in_interval[1], n-r, m_beta, m_eta, times);
	}
	cp = cp/B;

	return cp;
}

double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times)
{
	double value, prob;

	prob = condiprob(beta, eta, times);
	value = pbinom(upper, n, prob, 1, 0) - pbinom(lower, n, prob, 1, 0) + dbinom(lower, n, prob, 0);

	return value;
}

double find_binom_prob_type2( double lower, double upper, int n, double beta, double eta, double start, double duration)
{
	double prob, value;

	prob = condiprob_type2(beta, eta, start, duration);
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

int qualify_discrete_type2(double betab, double etab, double beta_mle, double eta_mle, double start, double duration){
	int flag1, flag2;
	double binomp;
	binomp = condiprob_type2(betab, etab, start, duration);
	flag1 = isnan(binomp);

	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	mub = log(etab);
	sigmab = 1/betab;

	betagpq = 1/(sigmah/sigmab*sigmah);
	etagpq = exp(muh + (muh-mub)/sigmab*sigmah);
	
	binomp = condiprob_type2(betab, etab, start, duration);

	flag2 = isnan(binomp);

	return (flag1 || flag2);
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
