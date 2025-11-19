#ifndef GAMMA_WEIGHT_FUNCTION
#define GAMMA_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "fntl.h"
#include "WeightFunction.h"

class GammaWeightFunction : public WeightFunction
{
public:
	GammaWeightFunction(double z, double alpha, double beta, double tol = 1e-6)
	: WeightFunction(), _z(z), _alpha(alpha), _beta(beta), _tol(tol)
	{
	}

	// Evaluate the function
	double eval(double x, bool log = false) const
	{
		double out;
		if (std::isinf(x) && x < 0) {
			out = -std::numeric_limits<double>::infinity();
		} else if (x >= _z) {
			out = -std::numeric_limits<double>::infinity();
		} else {
			out = (_alpha - 1) * std::log(_z-x) - (_z-x) / _beta + log(x < _z)-
				lgamma(_alpha) - _alpha * std::log(_beta);
		}

		return log ? out : std::exp(out);
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const
	{
		// Arguments for findroot
		fntl::findroot_args args;
		args.tol = _tol;
		
        double x1;
        double x2;
		double x_max = _z - _beta * (_alpha - 1);

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
				const fntl::dfd& pred = [&](double x) -> double {
					double out;
					if (std::isinf(x) && x < 0) {
						out = -std::numeric_limits<double>::infinity();
					} else if (x >= _z) {
						out = -std::numeric_limits<double>::infinity();
					} else {
						out = (_alpha - 1) * log(_z-x) - (_z-x) / _beta + log(x < _z) -
							lgamma(_alpha) - _alpha * log(_beta) - _log_a;
					}

					return out;
				}

				const auto& out1 = fntl::findroot_bisect(pred, x_lo, x_max, args);
				x1 = out1.root;

				const auto& out2 = fntl::findroot_bisect(pred, x_max, _z, args);
				x2 = out2.root;
			}
		}
		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		double x_max = _z - _beta * (_alpha - 1);
		return eval(x_max, true);
	}

private:
	double _z;
	double _alpha;
	double _beta;
	double _tol;
};

#endif
