Pulsatile
=========

[![Build Status](https://travis-ci.com/mmulvahill/pulsatile.svg?token=Vzy3B4WH2SvZ4ybN4Uzy&branch=master)](https://travis-ci.com/mmulvahill/pulsatile) [![codecov](https://codecov.io/gh/mmulvahill/pulsatile/branch/master/graph/badge.svg?token=WeMubsj4Is)](https://codecov.io/gh/mmulvahill/pulsatile)

An R package for analyzing time series of hormone concentrations using Bayesian deconvolution analysis.

Run the following code to install the development version:

``` r
# Need devtools package to install from GitHub
install.packages("devtools")
```

    ## Installing package into '/home/matt/R/x86_64-pc-linux-gnu-library/3.3'
    ## (as 'lib' is unspecified)

``` r
library(devtools)

# install from github
install_github("mmulvahill/pulsatile", auth_token)
```

    ## Using GitHub PAT from envvar GITHUB_PAT

    ## Downloading GitHub repo mmulvahill/pulsatile@master
    ## from URL https://api.github.com/repos/mmulvahill/pulsatile/zipball/master

    ## Installing pulsatile

    ## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
    ##   --quiet CMD INSTALL  \
    ##   '/tmp/RtmpVL9oCi/devtools1863576c4910/mmulvahill-pulsatile-ac7bc2cbdc4d01446a4a0b22cb7f9f46d525b4f8'  \
    ##   --library='/home/matt/R/x86_64-pc-linux-gnu-library/3.3'  \
    ##   --install-tests

    ## 

References
==========

\[R source (math.h and internal.h)\]<https://github.com/wch/r-source>

R-&gt;C-&gt;R
-------------

\[C interface in Advanced R\]<http://adv-r.had.co.nz/C-interface.html> \[R's RNG wrappers\]<https://svn.r-project.org/R/trunk/src/library/stats/R/distn.R> \[R's rmultinom\]<https://svn.r-project.org/R/trunk/src/nmath/rmultinom.c>
