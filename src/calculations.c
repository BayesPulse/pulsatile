//                                                                         
// calculations.c
//                                                                         
//   This routine takes a symmetric positive definite matrix and performs    
//   a cholesky decomposition.  It replaces the lower triangular part of     
//   A with G.                                                               
//                                                                           
//   The cholesky decomposition, decomposes A into A = GG'.  This version of 
//   the algorithm is an outer product (rank-1) update of the principal      
//   submatrices of A.                                                       
//                                                                           
//   See "Matrix Computations" by Golub and Van Loan.  2nd edition, Johns    
//   Hopkins University Press, 1989. p. 423ff                                
//                                                                           
//   Input  A - a SPD matrix                                                 
//          num_col - size of A                                              
//                                                                           
//   Returns 1 upon success, 0 if A not SPD                                  
// 

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "calculations.h"
//#include "deconvolution_main.h"



//-----------------------------------------------------------------------------
//
// cholesky_decomp
//
//-----------------------------------------------------------------------------
int cholesky_decomp(double **A, int num_col) {

  int k;
  int i;
  int j;
  long double tmp;

  /* this is the GAXPY update */
  for (j = 0; j < num_col; j++) {

    if (j > 0) {
      for (i = j; i < num_col; i++) {
        tmp = 0;
        for (k = 0; k <= j-1; k++) {
          tmp += A[i][k] * A[j][k];
        }
        A[i][j] -= tmp;
      }
      if (A[j][j] <= 0) {
        Rprintf("cholesky_decomp: error matrix A not SPD %d\n", j);
        Rprintf("j = %d %lf\n", j, A[j][j]);
        return 0;
      }
      tmp = sqrt(A[j][j]);
      for (i = j; i < num_col; i++) {
        A[i][j] /= tmp;
      }
    } else {
      tmp = sqrt(A[j][j]);
      for (i = j; i < num_col; i++)
        A[i][j] /= tmp;
    }

  }

  return 1;

}




//-----------------------------------------------------------------------------
// 
// **cholesky_invert
// 
//-----------------------------------------------------------------------------
double **cholesky_invert(int len, double **H) {

  // takes G from GG' = A and computes A inverse 
  int i; 
  int j; 
  int k;
  double temp;
  double **INV;

  INV = (double **)calloc(len, sizeof(double *));
  for (i = 0; i < len; i++) {
    INV[i] = (double *)calloc(len, sizeof(double));
  }
  for (i = 0; i < len; i++) {
    INV[i][i] = 1;
  }

  for (k = 0; k < len; k++) {

    INV[k][0] /= H[0][0];
    for (i = 1; i < len; i++) {
      temp = 0.0;
      for (j = 0; j < i; j++) {
        temp += H[i][j] * INV[k][j];
      }
      INV[k][i] = (INV[k][i] - temp)/H[i][i];
    }

    INV[k][len-1] /= H[len-1][len-1];
    for (i = len - 2; i >= 0; i--) {
      temp = 0.0;
      for (j = i + 1; j < len; j++) {
        temp += H[j][i] * INV[k][j];
      }
      INV[k][i] = (INV[k][i] - temp)/H[i][i];
    }

  }

  return INV;

}




//-----------------------------------------------------------------------------
// 
// rmvnorm
// 
//-----------------------------------------------------------------------------
int rmvnorm(double *result, 
            double **A, 
            int size_A, 
            double *mean, 
            //unsigned long *seed, 
            int flag) {

  int i;
  int j;
  double *runiv;

  for (i = 0; i < size_A; i++) {
    result[i] = 0;
  }

  j = 1;
  if (!flag) {
    j = cholesky_decomp(A, size_A);
  }

  if (j) {

    runiv = (double *)calloc(size_A, sizeof(double));
    for (i = 0; i < size_A; i++) {
      runiv[i] = Rf_rnorm(0, 1); //snorm(seed);
    }

    for (i=0;i<size_A;i++) {
      for (j=0;j<=i;j++) {
        result[i] += A[i][j]*runiv[j];
      }
      result[i] += mean[i];
    }

    free(runiv);
    return 1;

  } else {

    return 0;

  }

}



// Single random multinomial
int one_rmultinom(double *cumprobs, int n_probs) {

  //SEXP ans;
  //Rf_protect(ans = Rf_allocVector(INTSXP, n_probs));
  //probs = Rf_coerceVector(probs
 
  int i;
  int rtn = 0;
  int *ans;
  ans = (int *)calloc(n_probs, sizeof(int));
  double *probs;
  probs = (double *)calloc(n_probs, sizeof(double));

  for (i = 0; i < n_probs; i++) {
    if (i == 0) probs[i] = cumprobs[i];
    else probs[i] = cumprobs[i] - cumprobs[i-1];
    ans[i] = 0;
    //Rprintf("Probability for pulse %d = %f\n", i, probs[i]);
    //Rprintf("Ans for pulse %d = %d\n", i, ans[i]);
  }

  Rf_rmultinom(1, probs, n_probs, ans);

  for (i = 0; i < n_probs; i++) {
    if (ans[i] == 1) rtn = i;
    //Rprintf("What is ans[%d] after rmultinom?: %d\n", i, ans[i]);
  }
  //Rprintf("What is rtn after assignment?: %d\n", rtn);

  return(rtn);
}

////-----------------------------------------------------------------------------
////
//// rwishart
////
////-----------------------------------------------------------------------------
//int rwishart(double **result, double **S, int size_S, int df, unsigned long
//        *seed, int flag) {
//
//    int i,j,k;
//    double *x,*zero, **Stmp;
//    int rmvnorm(double *, double **, int, double *, unsigned long *, int);
//    
//    /*	if ((double)df < size_S) return 0; */
//
//    x    = (double *)calloc(size_S, sizeof(double));
//    zero = (double *)calloc(size_S, sizeof(double));
//    Stmp = (double **)calloc(size_S, sizeof(double *));
//
//    for (i = 0; i < 2; i++) {
//        Stmp[i] = (double *)calloc(size_S, sizeof(double));
//    }
//      
//    for (i = 0; i < size_S; i++) {
//        zero[i] = 0;
//    }
//    
//    for (i = 0; i < size_S; i++) {
//        for (j=0;j<size_S;j++) {
//            result[i][j] = 0;
//            Stmp[i][j] = S[i][j];  /* don't want to permanently alter the S
//                                      matrix if cholesky decomp not in main
//                                      program */
//        }
//    }
//
//    for (k = 0; k < df; k++) {
//        if (rmvnorm(x, Stmp, size_S, zero, seed, flag)) {
//            for (i = 0; i < size_S; i++) {
//                for (j = 0; j < size_S; j++) {
//                    result[i][j] += x[i] * x[j];
//                }
//            }
//        } else {
//            free(zero);
//            free(x);
//            for (i = 0; i < 2; i++) {
//                Stmp[i];
//            }
//            free(Stmp);
//            Rprintf("wishart asdf\n");
//            return 0;
//        }
//    }
//
//    free(zero);
//    free(x);
//
//    for (i = 0; i < 2; i++) {
//        Stmp[i];
//    }
//
//    free(Stmp);
//    return 1;
//
//}




//******************************************************************************
// END OF FILE
//******************************************************************************

