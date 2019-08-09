Pulsatile
=========

[![Build Status](https://travis-ci.org/BayesPulse/pulsatile.svg?branch=master)](https://travis-ci.org/BayesPulse/pulsatile)
[![codecov](https://codecov.io/gh/BayesPulse/pulsatile/branch/master/graph/badge.svg)](https://codecov.io/gh/BayesPulse/pulsatile)

An R package for analyzing time series of hormone concentrations using Bayesian deconvolution analysis.  This package is package is fully functional for fitting individual time series data.  It's primary purpose was for the development of the user interface in R for modeling this type of data and for implementing the R to C interface.  The backend algorithm, written in C, is functional, but not modularized or optimized.  That work is being done in the `libpulsatile` package.

Run the following code to install the development version:

``` r
# Need devtools package to install from GitHub
if (!require(devtools)) install.packages("devtools")

# install from github
install_github("bayespulse/pulsatile")
```

Programming resources
=====================

[R source (math.h and internal.h)](https://github.com/wch/r-source)
[C interface in Advanced R](http://adv-r.had.co.nz/C-interface.html)
[R's RNG wrappers](https://svn.r-project.org/R/trunk/src/library/stats/R/distn.R)
[R's rmultinom](https://svn.r-project.org/R/trunk/src/nmath/rmultinom.c)
