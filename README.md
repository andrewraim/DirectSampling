This package contains code to support the manuscript "Direct Sampling in
Bayesian Regression Models with Additive Disclosure Avoidance Noise" (link
pending). It may be extended beyond the manuscript as well.

To install the package from Github, first ensure that you have a recent
installation of R and the packages Rcpp, datastructures, R6, and devtools.
You must also have the ability to build C++ from within R; for example, Windows
can get [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac OS X
users can find similar tools
[here](https://cran.r-project.org/bin/macosx/tools/).

After installing the development tools, you should be able to run the following
following R script from within R and to print the vector `(2,3,4,5,6)`.

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

The folder `inst/manuscript` contains replication materials for the manuscript
which also demonstrate use of the package. See `inst/manuscript/README.md` for
a description of the contents.

