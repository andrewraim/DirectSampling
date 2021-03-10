---
title: DirectSampling Package
author: Andrew Raim
---

This package contains code to accompany the manuscript "Direct Sampling in
Bayesian Regression Models with Additive Disclosure Avoidance Noise" (link
pending).

To install the package from Github, first ensure that you have a recent
installation of R and the packages Rcpp, datastructures, R6, and devtools.
You must also have the ability to build C++ from within R; for example, Windows
users need Rtools <https://cran.r-project.org/bin/windows/Rtools/>. At a
minimum, the following R script should run and print the vector `(2,3,4,5,6)`.

```r
library(Rcpp)

Rcpp::sourceCpp(code = '
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector plus_one(const Rcpp::NumericVector& x) {
	return x + 1;
}
')

plus_one(1:5)
```

Now the DirectSampling package can be installed and built using devtools:
```
R> devtools::install_github("andrewraim/DirectSampling")
```

Package documentation is currently between minimal and nonexistent. However,
the folder `inst/manuscript` contains a README describing the replication
materials for the manuscript, which make use of the package.