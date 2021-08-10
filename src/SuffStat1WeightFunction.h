#ifndef SUFFSTAT1_WEIGHT_FUNCTION
#define SUFFSTAT1_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "WeightFunction.h"
#include "find_root.h"

/*
* This class is used in root-finding
*/
class SuffStat1Root_pred : public Functional1
{
public:
	SuffStat1Root_pred(double z, double sigma2, unsigned int n, unsigned int d, double log_a)
		: _z(z), _sigma2(sigma2), _n(n), _d(d), _log_a(log_a)
	{
	}

	double operator()(double x) const {
		double out;
		if (std::isinf(x) && x < 0) {
			out = -std::numeric_limits<double>::infinity();
		} else if (x >= _z) {
			out = -std::numeric_limits<double>::infinity();
		} else {
			out = 0.5*(_n-_d-2) * log(_z-x) - (_z-x) / (2*_sigma2) - _log_a;
		}
		// Rprintf("SuffStat1Root_pred(%g) = %g :: z = %g, sigma2 = %g, n = %d, d = %d, log_a = %g\n", x, out, _z, _sigma2, _n, _d, _log_a);
		return out;
	}

private:
	double _z;
	double _sigma2;
	unsigned int _n;
	unsigned int _d;
	double _log_a;
};

class SuffStat1WeightFunction : public WeightFunction
{
public:
	SuffStat1WeightFunction(double z, double sigma2, unsigned int n, unsigned int d, double tol = 1e-6)
		: WeightFunction(), _z(z), _sigma2(sigma2), _n(n), _d(d), _tol(tol)
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
			out = 0.5*(_n-_d-2) * log(_z-x) - (_z-x) / (2*_sigma2);
		}

		// Rprintf("SuffStat1WeightFunction.eval(%g) = %g :: z = %g, sigma2 = %g, n = %d, d = %d\n", x, out, _z, _sigma2, _n, _d);
		if (take_log) { return(out); } else { return(exp(out)); }
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const {
        double x1;
        double x2;
		double x_max = _z - _sigma2 * (_n-_d-2);

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
				SuffStat1Root_pred pred(_z, _sigma2, _n, _d, log_a);
				x1 = find_root(x_lo, x_max, pred, _tol);
				x2 = find_root(x_max, _z, pred, _tol);
			}
		}
		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		double x_max = _z - _sigma2 * (_n-_d-2);
		return eval(x_max, true);
	}
private:
	double _z;
	double _sigma2;
	unsigned int _n;
	unsigned int _d;
	double _tol;
};

#endif
