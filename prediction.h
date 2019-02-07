#define MAX_ARRAY 10000
#define N 10
#define B 5000

extern double censor;

double *single_continous_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper);
void single_binom_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper, double times);
double *simulator(int type, double Er, double Pt, double shape, double scale);
double *findmle(double data[], double weightArray[]);
double *generateWeights(int FRWB, int n);
double *bootSimulator(double Er, double Pt, double shape, double scale);
double gpqinterval(double betab[], double etab[], double alpha, double mbeta, double meta);
double pbinterval(double betab[], double etab[], double alpha);
double pbbinominterval(double betab[], double etab[], double alpha, double times, int n, int r);
double gpqbinominterval(double betab[], double etab[], double alpha, double times, double beta_mle, double eta_mle, int n, int r);
