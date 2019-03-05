#define MAX_ARRAY 6000
#define N 50
#define B 500

extern double censor;

double *single_continous_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper);
double * single_binom_iteration(int type, double Er, double Pt, double beta, double eta, int FRWB, double lower, double upper, double nextCen);
double *simulator(int type, double Er, double Pt, double shape, double scale);
double *findmle(double data[], double weightArray[]);
double *generateWeights(int FRWB, int r, int n);
double *bootSimulator(double Er, double Pt, double shape, double scale);
double *gpqinterval(double betab[], double etab[], double lower, double upper, double mbeta, double meta);
double *pbinterval(double betab[], double etab[], double lower, double upper);
double find_binom_prob(double lower, double upper, int n, double beta, double eta, double times);
double calibinominterval(double betab[], double etab[], double alpha, double times, double beta_mle, double eta_mle, int n, int r);

int *plug_in_binominterval(double beta, double eta, double lower, double upper, double times, int n, int r);
int *pbbinominterval(double betab[], double etab[], double lower, double upper, double times, int n, int r);
int *gpqbinominterval(double betab[], double etab[], double lower, double upper, double times, double beta_mle, double eta_mle, int n, int r);
int *fonsecabinominterval(double betab[], double etabt[], double lower, double upper, double times, double beta_mle, double eta_mle, int n, int r);
double single_binom_iteration_calibration_type1(double Er, double Pt, double beta, double eta, double lower, double upper, double nextCen);