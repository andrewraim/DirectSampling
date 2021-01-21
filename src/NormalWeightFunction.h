#ifndef NORMAL_WEIGHT_FUNCTION
#define NORMAL_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "WeightFunction.h"

class NormalWeightFunction : public WeightFunction
{
public:
	NormalWeightFunction(double z, double mu, double sigma2)
		: WeightFunction(), _z(z), _mu(mu), _sigma2(sigma2)
	{
	}
	// Evaluate the function
	double eval(double x, bool take_log = false) const {
		double out = -1/(2*_sigma2) * pow(_z - x - _mu, 2);
		if (take_log) { return(out); } else { return(exp(out)); }
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const {
        double x1;
        double x2;
        if (-2 * _sigma2 * log_a < 0 && log_a <= log_c()) {
			x1 = _z - _mu;
			x2 = _z - _mu;
        } else {
	        x1 = _z - (_mu + sqrt(-2 * _sigma2 * log_a));
	        x2 = _z - (_mu - sqrt(-2 * _sigma2 * log_a));
        }
		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		return 0;
	}
private:
	double _z;
	double _mu;
	double _sigma2;
};

#endif
