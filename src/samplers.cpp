#include "samplers.h"
#include "direct_sampler.h"
#include "direct_sampler_ar.h"
#include "LognormalWeightFunction.h"
#include "NormalWeightFunction.h"
#include "SuffStat1WeightFunction.h"
#include "SuffStat3WeightFunction.h"
#include "DGeomBaseDistribution.h"
#include "DscNormBaseDistribution.h"
#include "NormalBaseDistribution.h"
#include "LaplaceBaseDistribution.h"

Rcpp::NumericVector direct_sampler_lognormal_normal(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	LognormalWeightFunction w(z, mu, sigma2);
	NormalBaseDistribution g(tau);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}

Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DscNormBaseDistribution g(tau);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}

Rcpp::NumericVector direct_sampler_lognormal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	LognormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

Rcpp::NumericVector direct_sampler_lognormal_dgeom(unsigned int n, double z,
	double mu, double sigma2, double rho, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DGeomBaseDistribution g(rho);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}

Rcpp::NumericVector direct_sampler_normal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	NormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}

Rcpp::NumericVector direct_sampler_suffstat1_laplace(unsigned int n, double z,
	double sigma2, unsigned int m, unsigned int d, double lambda, double tol,
	unsigned int N, const std::string& fill_method, unsigned int max_rejections)
{
	SuffStat1WeightFunction w(z, sigma2, m, d);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}

Rcpp::NumericVector direct_sampler_suffstat3_laplace(unsigned int n, double z,
	double phi2, double sigma2, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections)
{
	SuffStat3WeightFunction w(z, phi2, sigma2);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method);
	}
}
