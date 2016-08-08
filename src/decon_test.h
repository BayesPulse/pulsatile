
#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

void decon_input(SEXP indata, SEXP model, SEXP iterations, SEXP thin);
double **convert_data(SEXP indata, int nrow); 
void deallocate_data(double **data, int nrow);
SEXP getListElement(SEXP list, const char *str);

