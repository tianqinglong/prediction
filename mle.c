#include <stdio.h>
#include <math.h>

#include "prediction.h"
#include "brent.h"

// make a local copy
static int r,n;
static double dataset[MAX_ARRAY], weight[MAX_ARRAY];

double geteta(double beta);
double func1(double beta);

double *findmle(double data[], double weightArray[])
{
	int i;
	static double MLEs[2];

	r = (int)data[0];
	n = (int)data[1];

	// Initialize the data and weights so that func1() can be used
	for(i=0;i<n;i++)
	{
		dataset[i] = data[i+2];
		weight[i] = weightArray[i];
	}

	// Find MLEs
	double machep = r8_epsilon();
	double t = machep;

	MLEs[0] = zero(0.001, 10000, machep, t, func1);
	MLEs[1] = geteta(MLEs[0]);

	return MLEs;
}

// This function gives the partial derivative that we want to solve.
// Use external variables "dataset" and "weight".
double func1(double beta)
{	
	int i;
	double value=0, eta;

// compute eta
	eta = geteta(beta);
	if(eta > 1000000000){
		eta = 1000000000;
	}

// compute the partial derivative
	for(i=0;i<n;i++)
	{
		value += ( i<r ) ? ( weight[i]/beta + weight[i]*log(dataset[i]) - weight[i]*log(eta) ) : 0;
		value -= weight[i] * pow(dataset[i]/eta,beta) * log(dataset[i]/eta);
	}

	return value;
}

//This function gives the mle of eta by providing the mle of beta
double geteta(double beta)
{
	int i;
	double eta;

	double nu=0, de=0;
	for(i=0;i<n;i++)
	{
		de += (i<r)?weight[i]:0;
		nu += weight[i]*pow(dataset[i],beta);
	}
	eta = pow(nu/de, 1/beta);

	return eta;
}
