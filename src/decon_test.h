
#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
//#include <stdlib.h>
//#include <math.h>
//#include <stdio.h>
//#include <time.h>

//SEXP decon_input(SEXP indata, SEXP model, SEXP iterations, SEXP thin);
double **convert_data(SEXP indata, int nrow); 
void deallocate_data(double **data, int nrow);
SEXP getListElement(SEXP list, const char *str);

typedef struct {
    double meanbh[2];    // prior mean on baseline and halflife
    double varbh[2];     // variance on prior of baseline and halflife
    double *fe_mean;     // array of prior mean for mean pulse mass and width
    double *fe_variance; // variance of prior for mean pulse mass
    double *re_sdmax;    // 'b' on uniform prior for random effects standard deviation
    double err_alpha;    // alpha: prior parameters for model error (inverse gamma)
    double err_beta;     // beta : ""
    double gamma;        // Gamma, repulsion parameter for Strauss prior
    double range;        // Range parameter for Strauss prior
} Priors;

typedef struct {
    double md[2];    // 0: baseline; 1: half-life
    double sigma;    // model error (variance)
    double lsigma;   // log of model error (may not be used)
    double theta[2]; // mean pulse mass and width on log scale
    double *re_sd;   // variance/covariance matrix for theta
    double decay;    // decay rate converted from above half-life
    double nprior;   // prior number of pulses
} Common_parms;

//#endif // DECONVOLUTION_MAIN_H
