#define ARRAY_MAX 100

extern int n, r;
extern double weiEta, weiBeta;
extern double cenArray[], realArray[], weightArray[];

void simuDataType2(double* cenVec, double* comVec);
double deriEquation(double weibullBeta);
double deriFRWB(double weibullBeta);
void getWeight(double *weight);
