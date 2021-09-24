#include "direct_sampler.h"
#include "Stepdown.h"

Rcpp::NumericVector direct_sampler(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N,
	const std::string& fill_method, double priority_weight)
{
	// Use our Stepdown approximation to draw from p(u)
	Stepdown step(w, g, tol, N, fill_method, priority_weight);
	const Rcpp::NumericVector& v = step.draw(n);
	const Rcpp::NumericVector& log_u = Rcpp::log(v);

	// Draw from g(x | u) for each value of log(u)
	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log_u(i));
		x(i) = g.r_truncated(endpoints.first, endpoints.second);
		Rcpp::checkUserInterrupt();
	}

	return x;
}
