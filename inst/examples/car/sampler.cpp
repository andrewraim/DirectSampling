// [[Rcpp::depends(DirectSampling,fntl)]]
#include "DirectSampling.h"
#include "BetaBaseDistribution.h"
#include "CARWeightFunction.h"

namespace ds = DirectSampling;

// [[Rcpp::export]]
Rcpp::List r_car_beta(unsigned int n, const Rcpp::NumericVector& lambda,
	double C, double x_lower, double x_upper, double a_rho, double b_rho,
	unsigned int N, unsigned int max_rejects, double tol = 1e-10)
{
	CARWeightFunction w(lambda, C, x_lower, x_upper, 1e-10);
	BetaBaseDistribution g(a_rho, b_rho);

	ds::stepdown_args args1;
	args1.N = N;
	args1.tol = tol;
	ds::Stepdown step(w, g, args1);

	ds::rejection_args args2;
	args2.max_rejects = max_rejects;
	const Rcpp::List& out = ds::rejection(step, n, args2);

	return out;
}

