
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

//#endif // DECONVOLUTION_MAIN_H
