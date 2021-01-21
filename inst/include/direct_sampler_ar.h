#ifndef DIRECT_SAMPLER_AR_H
#define DIRECT_SAMPLER_AR_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"

Rcpp::NumericVector direct_sampler_ar(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N_init,
	unsigned int max_rejections, const std::string& fill_method);
#endif
