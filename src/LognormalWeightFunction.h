#ifndef LOGNORMAL_WEIGHT_FUNCTION
#define LOGNORMAL_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "WeightFunction.h"

class LognormalWeightFunction : public WeightFunction
{
public:
	LognormalWeightFunction(double z, double mu, double sigma2)
		: WeightFunction(), _z(z), _mu(mu), _sigma2(sigma2)
	{
	}

	// Evaluate the function
	double eval(double x, bool take_log = false) const {
		double out = -INFINITY;
		if (x < _z) {
			out = -log(_z-x) - pow(log(_z-x) - _mu, 2) / (2*_sigma2);
		}
		if (take_log) { return(out); } else { return(exp(out)); }
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const {
        double x1;
        double x2;
        if (_sigma2 < 2*(_mu + log_a) && log_a <= log_c()) {
			x1 = _z - exp(_mu - _sigma2);
			x2 = _z - exp(_mu - _sigma2);
        } else {
			x1 = _z - exp((_mu - _sigma2) + sqrt(_sigma2 * (_sigma2 - 2*(_mu + log_a))));
			x2 = _z - exp((_mu - _sigma2) - sqrt(_sigma2 * (_sigma2 - 2*(_mu + log_a))));
        }

		// Edge case: If z is an integer and the larger root is numerically close to
		// z, make the root smaller by an epsilon. This helps to ensure that the
		// endpoint z is not part of the support of p(u).
		if (_z - floor(_z) < 1e-10 && _z - x2 < 1e-6) {
			x2 = _z - 1e-6;
		}
		x1 = std::min(x1, x2);

		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		return -(_mu - _sigma2 / 2);
	}
private:
	double _z;
	double _mu;
	double _sigma2;
};

#endif
