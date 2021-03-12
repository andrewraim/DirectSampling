#ifndef DSCNORM_H
#define DSCNORM_H

#include <Rcpp.h>

//' @name Discrete Normal
//' @export
// [[Rcpp::export]]
int hi_dscnorm(double sigma, double tol);

double d_dscnorm(int x, double sigma, double tol, bool take_log, bool normalize);
double p_dscnorm(int x, double sigma, double tol);
int r_dscnorm(double sigma);
int q_dscnorm(double q, double sigma, double tol);


//' Discrete Normal Distribution
//' 
//' Functions for the Discrete Normal distribution with mean zero.
//' 
//' @param n Number of draws to generate.
//' @param x A vector of points to evaluate.
//' @param q A vector of probabilities to evaluate.
//' @param sigma Standard deviation parameter.
//' @param take_log If \code{TRUE} return the log-density.
//' @param tol Tolerance to truncate to a finite support (see details).
//' @param normalize If \code{TRUE}, normalize the density (approximately).
//' 
//' @details
//' Here we assume Discrete Normal distribution with density
//' \deqn{
//' f(x) \propto \exp( -x^2 / (2 \sigma^2) ), \quad x \in
//' \{ \ldots, -1, 0, 1, \ldots \}.
//' }
//' 
//' The function \code{hi_dscnorm} returns an \eqn{x} so that \eqn{[-x, x]}
//' contains \code{1 - tol} mass of the distribution. This is used to truncate
//' the distribution to have finite support for the simple implementation here.
//'
//' @references
//' Clement L. Canonne, Gautam Kamath, Thomas Steinke (2021), The Discrete
//' Gaussian for Differential Privacy. <https://arxiv.org/abs/2004.00010>.
//'
//' @name Discrete Normal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector d_dscnorm(const Rcpp::NumericVector& x, double sigma,
	double tol, bool take_log, bool normalize);

//' @name Discrete Normal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector p_dscnorm(const Rcpp::NumericVector& x, double sigma, double tol);

//' @name Discrete Normal
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector r_dscnorm(unsigned int n, double sigma);

//' @name Discrete Normal
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector q_dscnorm(const Rcpp::NumericVector& q, double sigma, double tol);

#endif
