#ifndef FIND_INTERVAL_H
#define FIND_INTERVAL_H

#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
int find_interval(double x, const Rcpp::NumericVector& cutpoints);

#endif
