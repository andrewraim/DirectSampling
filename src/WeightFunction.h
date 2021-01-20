#ifndef WEIGHT_FUNCTION_H
#define WEIGHT_FUNCTION_H

#include <Rcpp.h>

class WeightFunction {
public:
	// Evaluate this function on x. If take_log, return the logarithm of the
	// value.
	virtual double eval(double x, bool take_log = false) const = 0;

	// Fine the roots of the equation log(w(x)) = log_a.
	virtual std::pair<double,double> roots(double log_a) const = 0;

	// Return the logarithm of maximum value of this function.
	virtual double log_c() const = 0;
};

#endif