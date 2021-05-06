#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

/*
* These didn't exist in Rcpp Sugar (as far as I could tell), so we'll make our
* own quick versions.
*/
Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);
Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease);
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x);

/*
* Compute quantiles of a discrete distribution with values 0, 1, ..., k-1 and
* associated *cumulative* probabilities cp(0), cp(1), ..., cp(k-1). Use a
* bisection search in case cp is a large vector . q and cp can be given on
* the log-scale or probability scale, but they are expected to be compatible.
*/
// [[Rcpp::export]]
unsigned int q_discrete(double q, const Rcpp::NumericVector& cp);

#endif
