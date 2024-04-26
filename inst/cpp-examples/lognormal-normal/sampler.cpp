// [[Rcpp::depends(DirectSampling)]]
#include "DirectSampling.h"
#include "../shared/LognormalWeightFunction.h"
#include "../shared/NormalBaseDistribution.h"

//' Specific Direct Sampler
//' 
//' Functions for a specific direct sampler.
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
//' Implementation of the following direct sampler:
//' \code{direct_sampler_lognormal_normal}: Lognormal weight function
//' \eqn{w(x \mid \mu, \sigma^2)} and Normal base distribution
//' \eqn{g(x \mid 0, \tau^2)}.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z,
	double mu, double sigma2, double tau, unsigned int N,
	unsigned int max_rejections, double tol = 1e-8)
{
	LognormalWeightFunction w(z, mu, sigma2);
	NormalBaseDistribution g(tau);
	const Rcpp::List& out = DirectSampling::direct_sampler_ar(n, w, g, tol, N,
		max_rejections, "small_rects", 0.5);
	return out;
}
