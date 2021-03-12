#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <Rcpp.h>

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
	double mu, double sigma2, double tau, double tol,
	unsigned int N, const std::string& fill_method);

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_dscnorm(unsigned int n, double z,
	double mu, double sigma2, double tau, double tol,
	unsigned int N, const std::string& fill_method);

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method);

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_lognormal_dgeom(unsigned int n, double z,
	double mu, double sigma2, double rho, double tol,
	unsigned int N, const std::string& fill_method);

//' @name Specific Direct Samplers
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector direct_sampler_normal_laplace(unsigned int n, double z,
	double mu, double sigma2, double lambda, double tol,
	unsigned int N, const std::string& fill_method);

#endif
