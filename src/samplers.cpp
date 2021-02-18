#include "samplers.h"
#include "direct_sampler.h"
#include "LognormalWeightFunction.h"
#include "NormalWeightFunction.h"
#include "DGeomBaseDistribution.h"
#include "DscNormBaseDistribution.h"
#include "NormalBaseDistribution.h"
#include "LaplaceBaseDistribution.h"

Rcpp::NumericVector direct_sampler_lognormal_normal(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol,
	unsigned int N, const std::string& fill_method)
{
	LognormalWeightFunction w(z, mu, sigma2);
	NormalBaseDistribution g(tau);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol,
	unsigned int N, const std::string& fill_method)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DscNormBaseDistribution g(tau);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

Rcpp::NumericVector direct_sampler_lognormal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method)
{
	LognormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

Rcpp::NumericVector direct_sampler_lognormal_dgeom(unsigned int n, double z,
	double mu, double sigma2, double rho, double tol,
	unsigned int N, const std::string& fill_method)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DGeomBaseDistribution g(rho);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

Rcpp::NumericVector direct_sampler_normal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method)
{
	NormalWeightFunction w(z, mu, sigma2);
	LaplaceBaseDistribution g(lambda);
	return direct_sampler(n, w, g, tol, N, fill_method);
}

