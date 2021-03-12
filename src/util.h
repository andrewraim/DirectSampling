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

#endif
