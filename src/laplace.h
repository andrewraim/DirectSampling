#ifndef LAPLACE_H
#define LAPLACE_H

#include <Rcpp.h>

double d_laplace(double x, double mu, double lambda, bool take_log);
double p_laplace(double x, double mu, double lambda);
double r_laplace(double mu, double lambda);
double q_laplace(double q, double mu, double lambda);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_laplace(const Rcpp::NumericVector& x, double mu, double lambda, bool take_log);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_laplace(const Rcpp::NumericVector& x, double mu, double lambda);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector r_laplace(unsigned int n, double mu, double lambda);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector q_laplace(const Rcpp::NumericVector& x, double mu, double lambda);

#endif
