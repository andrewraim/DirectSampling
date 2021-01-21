#ifndef DSCNORM_H
#define DSCNORM_H

#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
int hi_dscnorm(double sigma, double tol);

double d_dscnorm(int x, double sigma, double tol, bool take_log, bool normalize);
double p_dscnorm(int x, double sigma, double tol);
int r_dscnorm(double sigma);
int q_dscnorm(double q, double sigma, double tol);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_dscnorm(const Rcpp::NumericVector& x, double sigma,
	double tol = 1e-10, bool take_log = false, bool normalize = true);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_dscnorm(const Rcpp::NumericVector& x, double sigma, double tol = 1e-10);

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector r_dscnorm(unsigned int n, double sigma);

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector q_dscnorm(const Rcpp::NumericVector& q, double sigma, double tol = 1e-10);

#endif
