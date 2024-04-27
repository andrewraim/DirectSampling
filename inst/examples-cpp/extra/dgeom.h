#ifndef DGEOM_H
#define DGEOM_H

#include <Rcpp.h>

double d_dgeom(int x, double p, bool take_log = false);
double p_dgeom(double x, double p);
int r_dgeom(double p);
double q_dgeom(double q, double p);

//' Double Geometric Distribution
//' 
//' Functions for the Double Geometric distribution.
//' 
//' @param n Number of draws to generate.
//' @param x A vector of points to evaluate.
//' @param q A vector of probabilities to evaluate.
//' @param p probability parameter.
//' @param take_log If \code{TRUE} return the log-density.
//' 
//' @details
//' Here we assume Double Geometric distribution with density
//' \deqn{
//' f(x) = \frac{\rho}{2 - \rho} (1 - \rho)^{|x|}, \quad x \in
//' \{ \ldots, -1, 0, 1, \ldots \}.
//' }
//' 
//' @name Double Geometric
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_dgeom(const Rcpp::NumericVector& x, double p, bool take_log);

//' @name Double Geometric
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_dgeom(const Rcpp::NumericVector& x, double p);

//' @name Double Geometric
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector r_dgeom(unsigned int n, double p);

//' @name Double Geometric
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector q_dgeom(const Rcpp::NumericVector& q, double p);

#endif
