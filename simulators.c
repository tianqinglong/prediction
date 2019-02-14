#include <stdio.h>
#include <math.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "prediction.h"

double *simulator(int type, double Er, double Pt, double shape, double scale)
// This function can simulate both type 1 and type 2 data
{
	int i, r, n;
	double x, y, u;
	static double data[MAX_ARRAY];

	if(type == 1)
	// simulate type 1 data
	{
		r=0;
		n=Er/Pt;

		// make sure there are at least 2 failures
		while(r<2)
		{
			r=0;
			x=0;
			for(i=0;i<n;i++)
			{
				u = unif_rand();
				y = 1 - (1-x)*pow(1-u,1.0/(n-i));

				if(y >Pt)
				{
					break;
				}

				r++;
				data[i+2] = scale*pow(-log(1-y), 1/shape);
				x = y;
			}
		}
		
		for(i=r;i<n;i++)
		{
			data[i+2] = censor;
		}
		
		data[0] = r;
		data[1] = n;
	}
	else if(type == 2)
	// simulate type 2 data
	{
		r=Er;
		n=Er/Pt;

		x=0;
		for(i=0;i<r;i++)
		{
			u = unif_rand();
			y = 1 - (1-x)*pow(1-u,1.0/(n-i));

			data[i+2] = scale*pow(-log(1-y), 1/shape);
			x = y;
		}

		for(i=r;i<n;i++)
		{
			data[i+2] = data[r+1];
		}

		data[0] = r;
		data[1] = n;
	}
	return data;
}

double *bootSimulator(double Er, double Pt, double shape, double scale)
// This function will do type 1 censoring for the 
{
	int i;
	int r, n = Er/Pt;
	double x,y,u,temp;

	static double data[MAX_ARRAY];

	r=0;
	while(r<2)
	{
		r=0;x=0;
		for(i=0;i<n;i++)
		{
			u = unif_rand();
			y = 1 - (1-x)*pow(1-u,1.0/(n-i));

			temp = scale*pow(-log(1-y), 1/shape);
			if(temp > censor){
				break;
			}
			r++;
			x = y;
			data[i+2] = temp;
		}
	}

	for(i=r;i<n;i++)
	{
		data[i+2] = censor;
	}

	data[0] = r;
	data[1] = n;

	return data;
}

double *generateWeights(int FRWB, int n)
// This function gives the weight array.
// FRWB = 0: all weights are 1
{
	int i;

	static double weight[MAX_ARRAY];

	if(FRWB == 0)
	{
		i=0;
		while(i<n)
		{
			weight[i] = 1;
			i++;
		}
	}
	else if(FRWB == 1)
	{
		double V=0,Y[n];
		for(i=0;i<n;i++)
		{
			Y[i] = rexp(1);
			V += Y[i];
		}
		for(i=0;i<n;i++)
		{
			weight[i] = Y[i]/V;
		}
	}

	return weight;
}
