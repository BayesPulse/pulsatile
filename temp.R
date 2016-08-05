
library(pryr)
library(inline)

rinternals <- file.path(R.home("modules"), "Rinternals.c")
file.show(rinternals)

one_rmultinom <- cfunction(c(probs = "numeric", n_probs = "integer"), "

  
  int *ans;
  ans = (int)calloc(INTEGER(n_probs), (sizeof(int)));
  rmultinom(1, probs, n_probs, ans);

  return(ans);

  ")

probs   <- c(0.1, 0.3, 0.1, 0.3, 0.2)
n_probs <- length(probs)

one_rmultinom(probs, n_probs)

########################################
library(MCMCpack)
rinvgamma # = 1/rgamma()



