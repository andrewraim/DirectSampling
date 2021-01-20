#ifndef DGEOM_H
#define DGEOM_H

#include <Rcpp.h>

double d_dgeom(int x, double p, bool take_log = false);
double p_dgeom(double x, double p);
int r_dgeom(double p);
double q_dgeom(double q, double p);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_dgeom(const Rcpp::NumericVector& x, double p, bool take_log);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_dgeom(const Rcpp::NumericVector& x, double p);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector r_dgeom(unsigned int n, double p);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector q_dgeom(const Rcpp::NumericVector& x, double p);

#endif
