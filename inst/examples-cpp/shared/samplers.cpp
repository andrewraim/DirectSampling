// #include "samplers.h"
#include "direct_sampler.h"
#include "direct_sampler_ar.h"
#include "LognormalWeightFunction.h"
#include "NormalWeightFunction.h"
#include "GammaWeightFunction.h"
#include "DGeomBaseDistribution.h"
#include "DscNormBaseDistribution.h"
#include "NormalBaseDistribution.h"
#include "LaplaceBaseDistribution.h"
#include "TDistDFWeightFunction.h"
#include "UniformBaseDistribution.h"

//' Specific Direct Samplers
//' 
//' Functions for some specific direct samplers.
//' 
//' @param n Number of draws to generate.
//' @param z A vector of points to evaluate.
//' @param mu Mean parameter for weight function.
//' @param sigma2 standard deviation parameter for weight function.
//' @param tau Standard deviation parameter for base distribution.
//' @param lambda Scale parameter for base distribution.
//' @param rho Probability parameter for base distribution.
//' @param tol Tolerance for step function approximation in customized sampler.
//' @param N Number of knots to use in approximation for \eqn{p(u)}.
//' @param fill_method Knot selection method for customized direct sampler.
//' Can be either \code{equal_steps} or \code{small_rects}.
//' @param max_rejections Maximum number of allowed rejections in accept-reject.
//' If zero (the default), use step function to approximate \eqn{p(u)}, rather
//' than as an envelope; this will avoid any rejections but yield an approximate
//' sample.
//' 
//' @details
//' Implementations for the following direct samplers are provided:
//' \itemize{
//' \item \code{direct_sampler_lognormal_normal}: Lognormal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Normal base distribution
//' \eqn{g(x \mid 0, \tau^2)}.
//' 
//' \item \code{direct_sampler_lognormal_dscnorm}: Lognormal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Discrete Normal base distribution
//' \eqn{g(x \mid 0, \tau^2)}.
//' 
//' \item \code{direct_sampler_lognormal_laplace}: Lognormal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Laplace base distribution
//' \eqn{g(x \mid 0, \lambda)}.
//' 
//' \item \code{direct_sampler_lognormal_dgeom}: Lognormal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Double Geometric base distribution
//' \eqn{g(x \mid \rho)}.
//' 
//' \item \code{direct_sampler_normal_laplace}: Normal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Laplace base distribution
//' \eqn{g(x \mid 0, \lambda)}.
//' }
//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
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

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight, double dscnorm_tol)
{
	LognormalWeightFunction w(z, mu, sigma2);
	DscNormBaseDistribution g(tau, dscnorm_tol);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
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

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
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

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
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

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_gamma_laplace(unsigned int n, double z,
	double alpha, double beta, double lambda, double tol, unsigned int N,
	const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	GammaWeightFunction w(z, alpha, beta);
	LaplaceBaseDistribution g(lambda);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_tdist_unif(unsigned int n, 
	unsigned int m, double A, double B, double nu_min, double nu_max, double tol,
	unsigned int N, const std::string& fill_method, unsigned int max_rejections,
	double priority_weight)
{
	TDistDFWeightFunction w(m, A, B, nu_min, nu_max, 1e-10);
	UniformBaseDistribution g(nu_min, nu_max);
	if (max_rejections > 0) {
		return direct_sampler_ar(n, w, g, tol, N, max_rejections, fill_method, priority_weight);
	} else {
		return direct_sampler(n, w, g, tol, N, fill_method, priority_weight);
	}
}
