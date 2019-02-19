#include <stdio.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

void source(const char *name)
{
    SEXP e;

    PROTECT(e = lang2(install("source"), mkString(name)));
    R_tryEval(e, R_GlobalEnv, NULL);
    UNPROTECT(1);
}

double *R_getmle(double data[])
{
	int lendata;
    static double val[2];
	lendata = (int) data[1] + 3 + (int) data[0];

	SEXP arg;
	PROTECT(arg = allocVector(REALSXP, lendata));
    memcpy(REAL(arg), data, lendata * sizeof(double));

    // int i;

    // printf("%d\n", lendata);
    // for(i=0;i<lendata;i++)
    // {
    //     printf("%f ", REAL(arg)[i]);
    // }

    SEXP Weibullmle_call;
    PROTECT(Weibullmle_call = lang2(install("Weibull.mle"), arg));

    int errorOccurred;
    SEXP ret = R_tryEval(Weibullmle_call, R_GlobalEnv, &errorOccurred);

    val[0] = REAL(ret)[0];
    val[1] = REAL(ret)[1];

    UNPROTECT(2);

    return val;
}

double *findWeibullMLEs(double data[])
{
    static double mle[2];
    double *temp;

    source("getweibullmle.r");
    temp = R_getmle(data);

    mle[0] = temp[0];
    mle[1] = temp[1];

    return mle;
}
