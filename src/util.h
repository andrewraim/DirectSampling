#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

/*
* These didn't exist in Rcpp Sugar (as far as I could see), so we'll make our own
* quick versions.
*/
Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);
Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease);
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x);

// Currently assumes both roots are real
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector quadratic_roots(double A, double B, double C);

// Compute log(p[0] + ... + p[k-1]) given log(p[0]), ..., log(p[k-1]).
// Treat p[0] as the "normalizing" probability.
// Uses the method from
// https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation/subtraction
// for stable calculation on the log scale.
//' @export
// [[Rcpp::export]]
double logsumprobs(const Rcpp::NumericVector& logprob);

#endif
