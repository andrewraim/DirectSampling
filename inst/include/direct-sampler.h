#ifndef DIRECT_SAMPLER_H
#define DIRECT_SAMPLER_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"
#include "Stepdown.h"

namespace DirectSampling {

Rcpp::NumericVector direct_sampler(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N,
	const std::string& fill_method, double priority_weight = 0.5)
{
	// Use our Stepdown approximation to draw from p(u)
	Stepdown step(w, g, tol, N, fill_method, priority_weight);
	const Rcpp::NumericVector& log_u = step.draw(n, true);

	// Draw from g(x | u) for each value of log(u)
	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log_u(i));
		x(i) = g.r_truncated(endpoints.first, endpoints.second);
		Rcpp::checkUserInterrupt();
	}

	return x;
}

}

#endif
