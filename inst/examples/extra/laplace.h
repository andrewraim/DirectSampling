#ifndef LAPLACE_H
#define LAPLACE_H

#include <Rcpp.h>

double d_laplace(double x, double mu, double lambda, bool log);
double p_laplace(double x, double mu, double lambda);
double r_laplace(double mu, double lambda);
double q_laplace(double q, double mu, double lambda);

//' Laplace Distribution
//' 
//' Functions for the Laplace distribution.
//' 
//' @param n Number of draws to generate.
//' @param x A vector of points to evaluate.
//' @param q A vector of probabilities to evaluate.
//' @param mu Location parameter.
//' @param lambda Scale parameter.
//' @param log If \code{TRUE} return the log-density.
//' 
//' @details
//' Here we assume Laplace distribution with density
//' \deqn{
//' f(x) = \frac{1}{2\lambda} e^{-|x| / \lambda}.
//' }
//' 
//' @name Laplace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_laplace(const Rcpp::NumericVector& x, double mu, double lambda, bool log);

//' @name Laplace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_laplace(const Rcpp::NumericVector& x, double mu, double lambda);

//' @name Laplace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector r_laplace(unsigned int n, double mu, double lambda);

//' @name Laplace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector q_laplace(const Rcpp::NumericVector& q, double mu, double lambda);

#endif
