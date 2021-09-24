#ifndef DIRECT_SAMPLER_H
#define DIRECT_SAMPLER_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"

Rcpp::NumericVector direct_sampler(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N,
	const std::string& fill_method, double priority_weight = 0.5);

#endif
