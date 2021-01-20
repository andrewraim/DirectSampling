#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>
#include "functionals.h"

Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);
Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease);
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x);
double bisection(double x_lo, double x_hi,
	const Predicate& pred,
	const Functional2& mid,
	const Functional2& dist, double tol);
int find_interval(double x, const Rcpp::NumericVector& cutpoints);

Rcpp::NumericVector quadratic_roots(double A, double B, double C);
double logsumprobs(const Rcpp::NumericVector& logprob);

unsigned int rcateg(const Rcpp::NumericVector& p);
Rcpp::IntegerVector rcateg(unsigned int n, const Rcpp::NumericVector& p);
Rcpp::IntegerVector rcateg_mat(const Rcpp::NumericMatrix& P);
Rcpp::IntegerVector rcateg_logp(unsigned int n, const Rcpp::NumericVector& logprob);
Rcpp::IntegerVector rcateg_logp_mat(const Rcpp::NumericMatrix& logprob);

#endif
