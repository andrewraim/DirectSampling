#ifndef DIRECT_SAMPLING_BISECTION_H
#define DIRECT_SAMPLING_BISECTION_H

#include <Rcpp.h>
#include "typedefs.h"

namespace DirectSampling {

inline double bisection(double x_lo, double x_hi,
	const predicate_t& f, const midpoint_t& mid,
	const distance_t& dist, double tol)
{
	double x = mid(x_lo, x_hi);

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		bool ind = f(x);
		x_lo = ind*x_lo + (1-ind)*x;
		x_hi = ind*x + (1-ind)*x_hi;
		x = mid(x_lo, x_hi);
	}

	return x;
}

}

#endif
