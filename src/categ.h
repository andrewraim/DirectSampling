#ifndef CATEG_H
#define CATEG_H

#include <Rcpp.h>

unsigned int rcateg(const Rcpp::NumericVector& p);

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg(unsigned int n, const Rcpp::NumericVector& p);

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_mat(const Rcpp::NumericMatrix& P);

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_logp(unsigned int n, const Rcpp::NumericVector& logprob);
Rcpp::IntegerVector rcateg_logp_mat(const Rcpp::NumericMatrix& logprob);

#endif