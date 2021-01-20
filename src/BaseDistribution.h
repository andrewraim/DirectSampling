#ifndef BASE_DISTRIBUTION_H
#define BASE_DISTRIBUTION_H

#include <Rcpp.h>

class BaseDistribution {
public:
	// Evaluate the density of this distribution on x. If take_log, return the
	// log-density.
	virtual double density(double x, bool take_log = false) const = 0;
	
	// Return the probability of the open interval (x1,x2): P(x1 < X < x2)
	virtual double pr_interval(double x1, double x2) const = 0;

	// Return the p quantile of this distribution after truncating to interval
	// [x_min, x_max].
	virtual double q_truncated(double p, double x_min, double x_max) const = 0;

	// Take a draw from this distribution after truncating to interval
	// [x_min, x_max].
	double r_truncated(double x_min, double x_max) const {
		double u = R::runif(0, 1);
		return q_truncated(u, x_min, x_max);
	}
};

#endif
