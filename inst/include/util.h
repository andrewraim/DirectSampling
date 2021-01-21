#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);
Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease);
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x);

Rcpp::NumericVector quadratic_roots(double A, double B, double C);
double logsumprobs(const Rcpp::NumericVector& logprob);

#endif
