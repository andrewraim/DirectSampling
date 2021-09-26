#include "bisection.h"
#include <Rcpp.h>

double bisection(double x_lo, double x_hi,
	const Predicate& f,
	const Functional2& mid,
	const Functional2& dist, double tol)
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
