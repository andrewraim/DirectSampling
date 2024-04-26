#ifndef FIND_ROOT_H
#define FIND_ROOT_H

#include <Rcpp.h>
#include "bisection.h"

namespace DirectSampling {

double find_root(double x_lo, double x_hi, const Functional1& f, double tol)
{
	Midpoint mid;
	IntervalLength dist;
	double x = mid(x_lo, x_hi);

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		// If no sign change occurs between f(x) and f(x_lo), then take x_lo = x
		// else take x_hi = x
		bool ind = (f(x) >= 0) == (f(x_lo) >= 0);
		x_lo = ind*x + (1-ind)*x_lo;
		x_hi = ind*x_hi + (1-ind)*x;
		x = mid(x_lo, x_hi);
	}

	if (dist(x_lo, x_hi) > tol && (x_hi <= x_lo)) {
		Rcpp::stop("Numerical overflow in bisection. tol may be too small");
	}

	return x;

}

}

#endif
