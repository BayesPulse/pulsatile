///-----------------------------------------------------------------------------
///
/// FILE: randgen.h
///
/// DESCRIPTION: 
///   Function definitions for birthdeath_strauss.c
/// 
///-----------------------------------------------------------------------------

#ifndef RANDGEN_H
#define RANDGEN_H

double kiss(unsigned long *seed);

int runiform_n(int n, unsigned long *seed);

long rmultinomial(double *prob, long len, unsigned long *seed);

double runif_atob(unsigned long *seed, double a, double b);

double rexp(double beta, unsigned long *seed);

double sgamma(double a, unsigned long* seed);

double fsign( double num, double sign );

double snorm(unsigned long *seed);

double rnorm(double mean, double stdev, unsigned long *seed);

double rgamma(double alpha, double beta, unsigned long *seed);

double inverse_gamma(double alpha, double beta, unsigned long *seed);

double fact_ln(int k);

int rpois(double lambda, unsigned long *seed);

double *rdirichlet(double *alpha, int len, unsigned long *seed);

double rbeta(double alpha, double beta, unsigned long *seed);

#endif // RANDGEN_H

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

