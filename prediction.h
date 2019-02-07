#define MAX_ARRAY 10000
#define N 20
#define B 10000

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
double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times);
double calibinominterval(double betab[], double etab[], double alpha, double times, double beta_mle, double eta_mle, int n, int r);