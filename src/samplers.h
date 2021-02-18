#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <Rcpp.h>
#include "LognormalWeightFunction.h"
#include "NormalWeightFunction.h"
#include "DGeomBaseDistribution.h"
#include "DscNormBaseDistribution.h"
#include "LaplaceBaseDistribution.h"

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol,
	unsigned int N, const std::string& fill_method);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_dgeom(unsigned int n, double z,
	double mu, double sigma2, double rho, double tol,
	unsigned int N, const std::string& fill_method);

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_normal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method);

#endif
