// [[Rcpp::depends(DirectSampling)]]
#include "DirectSampling.h"
#include "BetaBaseDistribution.h"
#include "CARWeightFunction.h"

// [[Rcpp::export]]
Rcpp::List r_car_beta(unsigned int n, const Rcpp::NumericVector& lambda,
	double C, double x_lower, double x_upper, double a_rho, double b_rho,
	unsigned int N, unsigned int max_rejections, double tol = 1e-10)
{
	CARWeightFunction w(lambda, C, x_lower, x_upper, 1e-10);
	BetaBaseDistribution g(a_rho, b_rho);
	const Rcpp::List& out = DirectSampling::direct_sampler_ar(n, w, g, tol, N,
		max_rejections, "small_rects", 0.5);
	return out;
}

