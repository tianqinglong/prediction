#define ARRAY_MAX 200

extern int r,n;
extern double Pt, Er;
extern double weiEta, weiBeta;
extern double cenArray[], comArray[], weightArray[];

double deriEquation1(double weibullBeta);
void simuDataType1(double *cenVec, double *comVec);
void getWeight(double *weight);
double deriFRWB(double weibullBeta);