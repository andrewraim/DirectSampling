#ifndef FIND_INTERVAL_H
#define FIND_INTERVAL_H

#include <Rcpp.h>

//' Find Interval
//' 
//' @param x A number 
//' @param cutpoints A sorted vector
//' 
//' @details
//' Elements of cutpoints represent endpoints
//' of adjacent intervals \eqn{[c_0,c_1)}, \eqn{[c_1,c_2)} ...,
//' \eqn{[c_{k},c_{k+1})}. Return the index \code{i} such that \code{x} is in 
//' \eqn{[c_i,c_{i+1})}; or return \code{-1} if \eqn{x < c_0}
//' or \code{k+1} if \eqn{x > c_{k+1}}
//'
//' @export
// [[Rcpp::export]]
int find_interval(double x, const Rcpp::NumericVector& cutpoints);

#endif
