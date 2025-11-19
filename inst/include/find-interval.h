#ifndef DIRECT_SAMPLING_FIND_INTERVAL_H
#define DIRECT_SAMPLING_FIND_INTERVAL_H

#include <Rcpp.h>
#include "typedefs.h"
#include "bisection.h"

namespace DirectSampling {

//' Find Interval
//' 
//' @param x A number 
//' @param cutpoints A sorted vector
//' 
//' @details
//' Elements of cutpoints represent endpoints
//' of adjacent intervals \eqn{[c_0,c_1)}, \eqn{[c_1,c_2)} ...,
//' \eqn{[c_{k},c_{k+1})}. Return the index \code{i} such that \code{x} is in 
//' \eqn{[c_i,c_{i+1})}; or return \code{-1} if \eqn{x < c_0}
//' or \code{k+1} if \eqn{x > c_{k+1}}
//'
//' @export
// [[Rcpp::export]]
inline int find_interval(double x, const Rcpp::NumericVector& cutpoints)
{
	unsigned int N = cutpoints.size();
	unsigned int k = N-1;

	if (x < cutpoints(0)) {
		return -1;
	} else if (x > cutpoints(k)) {
		return N;
	}

	const predicate_t& pred = [&](double idx) -> bool {
		return cutpoints(int(idx)) >= x;
	};

	/* Midpoint with floor to get an integer */
	const midpoint_t& mid = [](double x, double y) -> double {
		return floor((y + x) / 2);
	};
	
	const distance_t& dist = [](double x, double y) -> double {
		return y - x;
	};
	
	double out = bisection(0, k, pred, mid, dist, 1);
	return int(out);
}

}

#endif
