#include "samplers.h"
#include "direct_sampler.h"
#include "direct_sampler_ar.h"
#include "LognormalWeightFunction.h"
#include "NormalWeightFunction.h"
#include "GammaWeightFunction.h"
#include "DGeomBaseDistribution.h"
#include "DscNormBaseDistribution.h"
#include "NormalBaseDistribution.h"
#include "LaplaceBaseDistribution.h"

Rcpp::NumericVector direct_sampler_lognormal_normal(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	LognormalWeightFunction w(z, mu, sigma2);
	NormalBaseDistribution g(tau);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DscNormBaseDistribution g(tau);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

Rcpp::NumericVector direct_sampler_lognormal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	LognormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else{
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

Rcpp::NumericVector direct_sampler_lognormal_dgeom(unsigned int n, double z,
	double mu, double sigma2, double rho, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DGeomBaseDistribution g(rho);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

Rcpp::NumericVector direct_sampler_normal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	NormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

Rcpp::NumericVector direct_sampler_gamma_laplace(unsigned int n, double z,
	double alpha, double beta, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	// Rprintf("direct_sampler_gamma_laplace: Checkpoint 0\n");
	GammaWeightFunction w(z, alpha, beta);
	// Rprintf("direct_sampler_gamma_laplace: Checkpoint 1\n");
	LaplaceBaseDistribution g(lambda);
	// Rprintf("direct_sampler_gamma_laplace: Checkpoint 2\n");
	if (max_rejections > 0) {
		// Rprintf("direct_sampler_gamma_laplace: Checkpoint 3a\n");
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		// Rprintf("direct_sampler_gamma_laplace: Checkpoint 3b\n");
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}
