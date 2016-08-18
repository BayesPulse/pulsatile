//
// FILE: cholesky.h
// DESCRIPTION: Function definitions for cholesky decomposition and random
//   multivariate normal
// 

#ifndef CALCULATIONS_H
#define CALCULATIONS_H

int cholesky_decomp(double **A, int num_col);

double **cholesky_invert(int len, double **H);

int rmvnorm(double *result, 
            double **A, 
            int size_A, 
            double *mean, 
            //unsigned long *seed, 
            int flag);

int one_rmultinom(double *probs, int n_probs);

//int rwishart(double **result, 
//             double **S, 
//             int size_S, 
//             int df, 
//             unsigned long *seed, 
//             int flag);

#endif 

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

