#ifndef SUFFSTAT3_WEIGHT_FUNCTION
#define SUFFSTAT3_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "WeightFunction.h"
#include "find_root.h"

/*
* This class is used in root-finding
*/
class SuffStat3Root_pred : public Functional1
{
public:
	SuffStat3Root_pred(double z, double phi2, double sigma2, double log_a)
		: _z(z), _phi2(phi2), _sigma2(sigma2), _log_a(log_a)
	{
	}

	double operator()(double x) const {
		double out;
		if (std::isinf(x) && x < 0) {
			out = -std::numeric_limits<double>::infinity();
		} else if (x >= _z) {
			out = -std::numeric_limits<double>::infinity();
		} else {
			out = 0.5 * log(_z - x) - _phi2 * (_z - x) / (2*_sigma2) - _log_a;
		}
		// Rprintf("SuffStat3Root_pred(%g) = %g :: z = %g, phi2 = %g, sigma2 = %g, log_a = %g\n", x, out, _z, _phi2, _sigma2, _log_a);
		return out;
	}

private:
	double _z;
	double _phi2;
	double _sigma2;
	double _log_a;
};

class SuffStat3WeightFunction : public WeightFunction
{
public:
	SuffStat3WeightFunction(double z, double phi2, double sigma2, double tol = 1e-6)
		: WeightFunction(), _z(z), _phi2(phi2), _sigma2(sigma2), _tol(tol)
	{
	}
	// Evaluate the function
	double eval(double x, bool take_log = false) const {
		double out;
		if (std::isinf(x) && x < 0) {
			out = -std::numeric_limits<double>::infinity();
		} else if (x >= _z) {
			out = -std::numeric_limits<double>::infinity();
		} else {
			out = 0.5 * log(_z - x) - _phi2 * (_z - x) / (2*_sigma2);
		}

		if (take_log) { return(out); } else { return(exp(out)); }
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const {
        double x1;
        double x2;
        double x_max = _z - _sigma2 / _phi2;

   		if (log_a >= log_c()) {
   			x1 = x_max;
   			x2 = x_max;
   		} else if (log_a < 0 && std::isinf(log_a)) {
   			x1 = -std::numeric_limits<double>::infinity();
   			x2 = _z;
		} else {
			// We probably need to use numerical root-finding, unless we're in
			// a special situation. First find an x_lo < x_max where
			// log[w(x_lo)] < log(a)
			double b = 0;
			double x_lo = x_max - pow(2, b);
			while (eval(x_lo, true) >= log_a) {
				b++;
				x_lo = x_max - pow(2, b);
			}

			if (std::isinf(x_lo)) {
				// If x_lo is too small to capture numerically, just return the whole support
	   			x1 = -std::numeric_limits<double>::infinity();
				x2 = _z;
			} else {
				// Use numerical root-finding
				SuffStat3Root_pred pred(_z, _phi2, _sigma2, log_a);
				x1 = find_root(x_lo, x_max, pred, _tol);
				x2 = find_root(x_max, _z, pred, _tol);
			}
		}
		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		double x_max = _z - _sigma2 / _phi2;
		return eval(x_max, true);
	}
private:
	double _z;
	double _phi2;
	double _sigma2;
	double _tol;
};

#endif
