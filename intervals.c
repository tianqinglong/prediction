#include <stdio.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <math.h>

#include "prediction.h"
#include "brent.h"

static double mlebeta, mleeta, level, betabt[B], etabt[B];

double funcpb(double x);
double funcgpq(double x);
double condiprob(double shape, double scale, double times);
double condiprob_type2( double shape, double scale, double censoring, double duration);
double calibrate(double level_cali);

/* The Continuous Cases */

// gpq method
double *gpqinterval(double betab[], double etab[], double lower, double upper, double mbeta, double meta)
{
	int i;
	double interval;
	static double gpqinterval[2];
	// make a static global copy
	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}

	mlebeta = mbeta;
	mleeta = meta;

	double machep = r8_epsilon();
	double t = machep;

	level = lower;
	interval = zero(0, 10000000, machep, t, funcgpq);
	gpqinterval[0] = interval;

	level = upper;
	interval = zero(0, 10000000, machep, t, funcgpq);
	gpqinterval[1] = interval;

	return gpqinterval;
}	

// percentile bootstrap
double *pbinterval(double betab[], double etab[], double lower, double upper)
{
	int i;
	static double pbinterval[2];

	double machep = r8_epsilon();
	double t = machep;

	// make a static global copy
	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}

	level = lower;
	pbinterval[0] = zero(0, 10000000, machep, t, funcpb);

	level = upper;
	pbinterval[1] = zero(0, 10000000, machep, t, funcpb);

	return pbinterval;
}


/* The Discrete Cases */

///////////////////////////
//Plug-in Method///////////
///////////////////////////

int *plug_in_binominterval(double beta, double eta, double lower, double upper, double times, int n, int r)
{
	static int interval[2];
	int temp;
	double binom_prob;
	binom_prob = condiprob(beta, eta, times);

	temp = qbinom(lower, n-r, binom_prob, 1, 0);
	interval[0] = ( temp == 0 ) ? 0 : temp - 1;

	temp = qbinom(upper, n-r, binom_prob, 1, 0);
	interval[1] = temp;

	return interval;
}

int *plug_in_binominterval_type2(double beta, double eta, double lower, double upper, double start, double duration, int n, int r)
{
	static int interval[2];
	int temp;
	double binom_prob;
	binom_prob = condiprob_type2( beta, eta, start, duration );

	// printf("%f\n", binom_prob);

	temp = qbinom(lower, n-r, binom_prob, 1, 0);
	interval[0] = ( temp == 0 ) ? 0 : temp - 1;

	temp = qbinom(upper, n-r, binom_prob, 1, 0);
	interval[1] = temp;

	return interval;
}

///////////////////////////
//Percentile Bootstrap/////
///////////////////////////

int *pbbinominterval(double betab[], double etab[], double lower, double upper, double times, int n, int r)
// percentile bootstrap for within sample discrete prediction
{
	int i, j;
	int lb = -10, ub = -10;
	double binomp, cdf, cdf_pre = 0;

	static int interval[2];

	for(i=0;i<=(n-r);i++)
	{
		cdf = 0;
		for(j=0;j<B;j++)
		{
			binomp = condiprob(betab[j], etab[j], times);
			cdf += pbinom(i,n - r,binomp,1,0);
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower)
		{
			lb = (i == 0) ? 0: (i-1);
		}

		if(cdf >= upper)
		{
			ub = i;
			break;
		}

		cdf_pre = cdf;
	}

	interval[0] = lb;
	interval[1] = ub;

	return interval;
}

int *pbbinominterval_type2(double betab[], double etab[], double censor_time[], double duration, double lower, double upper, int n, int r)
// percentile bootstrap for within sample discrete prediction
{
	int i, j;
	int lb = -10, ub = -10;
	double binomp, cdf, cdf_pre = 0;

	static int interval[2];

	for(i=0;i<=(n-r);i++)
	{
		cdf = 0;
		for(j=0;j<B;j++)
		{
			binomp = condiprob_type2(betab[j], etab[j], censor_time[j], duration);
			cdf += pbinom(i,n - r,binomp,1,0);
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower)
		{
			lb = (i == 0) ? 0: (i-1);
		}

		if(cdf >= upper)
		{
			ub = i;
			break;
		}

		cdf_pre = cdf;
	}

	interval[0] = lb;
	interval[1] = ub;

	return interval;
}
///////////////////////////
//GPQ Binom Interval///////
///////////////////////////

int *gpqbinominterval(double betab[], double etab[], double lower, double upper, double times, double beta_mle, double eta_mle, int n, int r)
{
	int i,j;
	double binomp, cdf, cdf_pre;
	int lb =-10, ub=-10;
	static int gpq_binom_interval[2];

	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	cdf_pre = 0;
	for(i=0;i<=(n-r);i++)
	{
		cdf=0;
		for(j=0;j<B;j++)
		{
			mub = log(etab[j]);
			sigmab = 1/betab[j];

			betagpq = 1/(sigmah/sigmab*sigmah);
			etagpq = exp(muh + (muh-mub)/sigmab*sigmah);

			binomp = condiprob(betagpq, etagpq, times);
			cdf += pbinom(i, n-r, binomp, 1, 0);
			if(isnan(binomp)){
				printf("%f %f %f %f %f %f\n", betagpq, etagpq, times, censor, betab[j], etab[j]);
				printf("%f %f\n", beta_mle, eta_mle);
			}
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower)
		{
			lb = ( i == 0 ) ? 0 : i-1;
		}

		if(cdf >= upper){
			ub = i;
			break;
		}

		cdf_pre = cdf;
	}

	if(lb == -10){
		printf("%f\n", cdf);
		getchar();
	}

	gpq_binom_interval[0] = lb;
	gpq_binom_interval[1] = ub;

	return gpq_binom_interval;
}

int *gpqbinominterval_type2(double betab[], double etab[], double censor_time[], double lower, double upper, double duration, double beta_mle, double eta_mle, int n, int r)
{
	int i,j;
	double binomp, cdf, cdf_pre;
	int lb =-10, ub=-10;
	static int gpq_binom_interval[2];

	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	cdf_pre = 0;
	for(i=0;i<=(n-r);i++)
	{
		cdf=0;
		for(j=0;j<B;j++)
		{
			mub = log(etab[j]);
			sigmab = 1/betab[j];

			betagpq = 1/(sigmah/sigmab*sigmah);
			etagpq = exp(muh + (muh-mub)/sigmab*sigmah);

			binomp = condiprob_type2(betagpq, etagpq, censor_time[j], duration);
			cdf += pbinom(i, n-r, binomp, 1, 0);
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower)
		{
			lb = ( i == 0 ) ? 0 : i-1;
		}

		if(cdf >= upper){
			ub = i;
			break;
		}

		cdf_pre = cdf;
	}

	gpq_binom_interval[0] = lb;
	gpq_binom_interval[1] = ub;

	return gpq_binom_interval;
}
///////////////////////////
//Fonseca Binom Interval///
///////////////////////////

int *fonsecabinominterval(double betab[], double etabt[], double lower, double upper, double times, double beta_mle, double eta_mle, int n, int r)
{
	int i, j, lb = -10, ub = -10;
	double cdf, pbm, pbb, cdf_pre=0;
	static int fonseca_interval[2];

	pbm = condiprob( beta_mle, eta_mle, times );
	
	for( i=0 ;i <= (n-r); i++ )
	{	
		cdf = 0;
		for(j=0;j<B;j++)
		{
			pbb = condiprob(betab[j], etabt[j], times);
			cdf += pbinom(qbinom(pbinom(i,n-r,pbm,1,0),n-r,pbb,1,0),n-r,pbm,1,0);
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower){
			lb = (i ==0) ? 0 : i-1;
		}

		if(cdf >= upper){
			ub = i;
			break;
		}
		cdf_pre = cdf;
	}

	fonseca_interval[0] = lb;
	fonseca_interval[1] = ub;

	return fonseca_interval;
}

int *fonsecabinominterval_type2(double betab[], double etabt[], double censor_time[], double lower, double upper, double duration, double beta_mle, double eta_mle, int n, int r, double original_censor)
{
	int i, j, lb = -10, ub = -10;
	double cdf, pbm, pbb, cdf_pre=0;
	static int fonseca_interval[2];

	pbm = condiprob_type2( beta_mle, eta_mle, original_censor, duration);
	
	for( i=0 ;i <= (n-r); i++ )
	{	
		cdf = 0;
		for(j=0;j<B;j++)
		{
			pbb = condiprob_type2(betab[j], etabt[j], censor_time[j], duration);
			cdf += pbinom(qbinom(pbinom(i,n-r,pbm,1,0),n-r,pbb,1,0),n-r,pbm,1,0);
		}
		cdf = cdf/B;

		if(cdf > lower && cdf_pre <= lower){
			lb = (i ==0) ? 0 : i-1;
		}

		if(cdf >= upper){
			ub = i;
			break;
		}
		cdf_pre = cdf;
	}

	fonseca_interval[0] = lb;
	fonseca_interval[1] = ub;

	return fonseca_interval;
}

// int *calibinominterval(double lower, double upper, double times, double beta_mle, double eta_mle, int n, int r)
// {

// }

///////////////////////////////
/* Below are local functions */
///////////////////////////////

// double (double beta, double eta)
// {

// 	return value;
// }

double condiprob(double shape, double scale, double times)
{
	double value;

	value = (pweibull(times*censor, shape, scale, 1, 0) - pweibull(censor, shape, scale, 1, 0))/(1-pweibull(censor, shape, scale, 1, 0));

	return value;
}

double condiprob_type2( double shape, double scale, double censoring, double duration)
{
	double value;

	double t2;
	double t1;

	t1 = censoring;
	t2 = censoring + duration;

	value = ( pweibull( t2, shape, scale, 1, 0) - pweibull(t1, shape, scale, 1, 0) )/( 1 - pweibull(t1, shape, scale, 1, 0) );

	return value;
}

double funcgpq(double x)
{
	int i;
	double cdf=0;
	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(mleeta);
	sigmah = 1/mlebeta;

	for(i=0;i<B;i++)
	{
		mub = log(etabt[i]);
		sigmab = 1/betabt[i];

		betagpq = 1/(sigmah/sigmab*sigmah);
		etagpq = exp(muh + (muh-mub)/sigmab*sigmah);

		cdf += pweibull(x, betagpq, etagpq, 1, 0);
	}
	cdf = cdf/B - level;

	return cdf;
}

double funcpb(double x)
{
	int i;
	double cdf=0;
	for(i=0;i<B;i++){
		cdf += pweibull(x, betabt[i], etabt[i], 1, 0);
	}
	cdf = cdf/B - level;

	return cdf;
}