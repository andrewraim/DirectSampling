#ifndef DIRECT_SAMPLING_DRAW_H
#define DIRECT_SAMPLING_DRAW_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"
#include "Stepdown.h"

namespace DirectSampling {

/*
* Draw from target using the step function. This will yield draws which are
* distributed approximately like the target.
* 
* Arguments
* - step: a given step function.
* - n: number of draws to generate.
*/
inline Rcpp::NumericVector draw(const Stepdown& step, unsigned int n)
{
	// Use our Stepdown approximation to draw from p(u)
	const Rcpp::NumericVector& log_u = step.draw(n, true);

	const WeightFunction& w = step.get_weight();
	const BaseDistribution& g = step.get_base();

	// Draw from g(x | u) for each value of log(u)
	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log_u(i));
		x(i) = g.r_truncated(endpoints.first, endpoints.second);

		if (i % 10000 == 0) {
		 	Rcpp::checkUserInterrupt();
		}
	}

	return x;
}

}

#endif
