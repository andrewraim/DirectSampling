#ifndef DSCNORM_H
#define DSCNORM_H

#include <Rcpp.h>

double d_dscnorm_unnormalized(int x, double sigma);
int r_dscnorm(double sigma);

// [[Rcpp::export]]
Rcpp::NumericVector d_dscnorm_unnormalized(const Rcpp::IntegerVector& x, double sigma);

// [[Rcpp::export]]
Rcpp::IntegerVector r_dscnorm(unsigned int n, double sigma);

#endif
