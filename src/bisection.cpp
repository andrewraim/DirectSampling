#include "bisection.h"

double bisection(double x_lo, double x_hi,
	const Predicate& pred,
	const Functional2& mid,
	const Functional2& dist, double tol)
{
	double x = mid(x_lo, x_hi);

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		bool ind = pred(x);
		x_lo = ind*x_lo + (1-ind)*x;
		x_hi = ind*x + (1-ind)*x_hi;
		x = mid(x_lo, x_hi);
	}

	if (dist(x_lo, x_hi) > tol && (x <= x_lo || x >= x_hi)) {
		Rcpp::stop("Numerical overflow in bisection. tol may be too small");
	}

	return x;
}
