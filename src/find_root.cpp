#include "find_root.h"

double find_root(double x_lo, double x_hi, const Functional1& f, double tol)
{
	Midpoint mid;
	IntervalLength dist;
	double x = mid(x_lo, x_hi);

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		// If no sign change occurs between f(x) and f(x_lo), then take x_lo = x
		// else take x_hi = x
		// Rprintf("Before: x_lo = %g, x_hi = %g, x = %g, f(x_lo) = %g, f(x_hi) = %g, f(x) = %g, tol = %g\n",
		// 	x_lo, x_hi, x, f(x_lo), f(x_hi), f(x), tol);
		bool ind = (f(x) >= 0) == (f(x_lo) >= 0);
		x_lo = ind*x + (1-ind)*x_lo;
		x_hi = ind*x_hi + (1-ind)*x;
		x = mid(x_lo, x_hi);
		// Rprintf("After: ind = %d, x_lo = %g, x_hi = %g, x = %g, f(x_lo) = %g, f(x_hi) = %g, f(x) = %g, tol = %g\n",
		// 	ind, x_lo, x_hi, x, f(x_lo), f(x_hi), f(x), tol);
	}

	if (dist(x_lo, x_hi) > tol && (x_hi <= x_lo)) {
		// Rprintf("x_lo = %g, x_hi = %g, x = %g, f(x_lo) = %g, f(x_hi) = %g, f(x) = %g, tol = %g\n",
		//	x_lo, x_hi, x, f(x_lo), f(x_hi), f(x), tol);
		Rcpp::stop("Numerical overflow in bisection. tol may be too small");
	}

	// Rprintf("Found root f(%g) = %g with x_lo = %g and x_hi = %g, dist(x_lo, x_hi) = %g\n", x, f(x), x_lo, x_hi, dist(x_lo, x_hi));
	// Rcpp::stop("PAUSE!");

	return x;

}
