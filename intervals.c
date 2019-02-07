#include <stdio.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <math.h>

#include "prediction.h"
#include "brent.h"

static double mlebeta, mleeta, level, betabt[B], etabt[B], local_time, local_n, local_r;

double funcpb(double x);
double funcgpq(double x);
double condiprob(double shape, double scale, double times);

double gpqinterval(double betab[], double etab[], double alpha, double mbeta, double meta)
// gpq method
{
	int i;
	double interval;
	// make a static global copy
	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}
	level = alpha;

	mlebeta = mbeta;
	mleeta = meta;

	double machep = r8_epsilon();
	double t = machep;

	interval = zero(0, 10000000, machep, t, funcgpq);

	return interval;
}

double pbinterval(double betab[], double etab[], double alpha)
// percentile bootstrap
{
	int i;
	double interval;

	double machep = r8_epsilon();
	double t = machep;

	// make a static global copy
	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}
	level = alpha;

	interval = zero(0.00001, 10000000, machep, t, funcpb);

	return interval;
}

double pbbinominterval(double betab[], double etab[], double alpha, double times, int n, int r)
// percentile bootstrap for within sample discrete prediction
{
	int i,j,bound;
	double binomp, cdf;
	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}

	level = alpha;
	local_time = times;

	local_r = r;
	local_n = n;

	for(i=0;i<=(n-r);i++)
	{
		cdf = 0;
		for(j=0;j<B;j++)
		{
			binomp = condiprob(betabt[j], etabt[j], local_time);
			cdf += pbinom(i,local_n-local_r,binomp,1,0);
		}
		cdf = cdf/B;

		if(cdf > level)
		{
			bound = i;
			break;
		}
	}

	if(level < 0.5)
	{
		bound = (bound==0)?0:(bound-1);
	}

	return bound;
}

double gpqbinominterval(double betab[], double etab[], double alpha, double times, double beta_mle, double eta_mle, int n, int r)
{
	int i,j,bound;
	double binomp, cdf;

	double sigmab, mub, sigmah, muh;
	double betagpq, etagpq;

	muh = log(eta_mle);
	sigmah = 1/beta_mle;

	for(i=0;i<B;i++)
	{
		betabt[i] = betab[i];
		etabt[i] = etab[i];
	}

	level = alpha;
	local_time = times;

	local_n = n;
	local_r = r;

	for(i=0;i<=(n-r);i++)
	{
		cdf=0;
		for(j=0;j<B;j++)
		{
			mub = log(etab[j]);
			sigmab = 1/betab[j];

			betagpq = 1/(sigmah/sigmab*sigmah);
			etagpq = exp(muh + (muh-mub)/sigmab*sigmah);

			binomp = condiprob(betagpq, etagpq, local_time);
			cdf += pbinom(i, local_n-local_r, binomp, 1, 0);
		}
		cdf = cdf/B;

		if(cdf > level)
		{
			bound = i;
			break;
		}
	}
	if(level < 0.5){
		bound = (bound == 0)?0:(bound-1);
	}

	return bound;
}

///////////////////////////////
/* Below are local functions */
///////////////////////////////

double condiprob(double shape, double scale, double times)
{
	double value;

	value = (pweibull(times*censor, shape, scale, 1, 0) - pweibull(censor, shape, scale, 1, 0))/(1-pweibull(censor, shape, scale, 1, 0));

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

// void debugger(double betab[], double etab[])
// {
// 	int i;
// 	//double x=0.1;

// 	// make a static global copy
// 	for(i=0;i<B;i++)
// 	{
// 		betabt[i] = betab[i];
// 		etabt[i] = etab[i];
// 	}
// 	// printf("x\tcdf(x)\n");
// 	// for(i=0;i<100;i++)
// 	// 	printf("%f %f\n", by*i, funcpb(by*i)+level);
// 	// for(i=0;i<B;i++)
// 	// {
// 	// 	printf("%f\n", pweibull(x,betabt[i],etabt[i],1,0));
// 	// }
// 	printf("%f %f\n", betabt[420],etabt[420]);
// }

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

