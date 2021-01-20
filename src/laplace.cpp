#include "laplace.h"

double d_laplace(double x, double mu, double lambda, bool take_log)
{
	double log_fx = -log(2) - log(lambda) -fabs(x - mu) / lambda;
	if (take_log) {
		return log_fx;
	} else {
		return exp(log_fx);
	}
}

double p_laplace(double x, double mu, double lambda)
{
	if (x <= mu) {
		return 0.5 * exp((x - mu) / lambda);
	} else {
		return 1 - 0.5 * exp(-(x - mu) / lambda);
	}
}

double r_laplace(double mu, double lambda)
{
	return R::rexp(lambda) - R::rexp(lambda) + mu;
}

double q_laplace(double q, double mu, double lambda)
{
	if (q <= 0.5) {
		return mu + lambda * log(2*q);
	} else {
		return mu - lambda * log(2 - 2*q);
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector d_laplace(const Rcpp::NumericVector& x, double mu, double lambda, bool take_log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_laplace(x(i), mu, lambda, take_log);
	}

	return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_laplace(const Rcpp::NumericVector& x, double mu, double lambda)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_laplace(x(i), mu, lambda);
	}

	return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector r_laplace(unsigned int n, double mu, double lambda)
{
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_laplace(mu, lambda);
	}

	return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector q_laplace(const Rcpp::NumericVector& x, double mu, double lambda)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_laplace(x(i), mu, lambda);
	}

	return out;
}
