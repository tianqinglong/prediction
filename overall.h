#define ARRAY_MAX 10000
#define B 5000

extern int n, r;
extern double Pt, Er;
extern double weiEta, weiBeta, alpha;
extern double cenArray[], realArray[], weightArray[];
extern double cenBootSample[], comBootSample[];
extern double betaBoot[], etaBoot[];
extern double etaMLE, betaMLE;

void simuDataType1(double *cenVec, double *comVec);
void simuDataType2(double* cenVec, double* comVec, double shape, double scale);
double deriEquation(double weibullBeta);
double deriFRWB(double weibullBeta);
void getWeight(double *weight);
double getIntervalsGPQ(double quantile);
double getIntervalsPercentileBoot(double quantile);
double getIntervalsFonseca(double quantile);
double deriEquationParaBoot(double shape);
