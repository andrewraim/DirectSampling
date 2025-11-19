// [[Rcpp::depends(DirectSampling)]]
#ifndef TDISTDF_WEIGHT_FUNCTION
#define TDISTDF_WEIGHT_FUNCTION

#include "DirectSampling.h"

enum TDistDFShape { null, mode, increasing, decreasing };

class TDistDFWeightFunction : public DirectSampling::WeightFunction
{
public:
	TDistDFWeightFunction(double n, double A, double nu_min, double nu_max,
		double tol)
	: DirectSampling::WeightFunction(), _n(n), _A(A), _nu_min(nu_min),
	  _nu_max(nu_max), _tol(tol), _grad_nu_min(0), _grad_nu_max(0),
	  _x_min(nu_min), _x_max(nu_max), _log_min(R_NegInf),
	  _log_max(R_NegInf), _shape(null)
	{
		init();
	}

	void init()
	{
		const DirectSampling::dfd_t& grad = [&](double x) -> double {
			return _n/2 * log(x/2) - _n/2 * R::digamma(x/2) + _n/2 - _A;
		};

		_grad_nu_min = grad(_nu_min);
		_grad_nu_max = grad(_nu_max);

		if (_n/2 - _A < 0) {
			// The condition A > 0 means the log w(x) is unimodal and the
			// mode is somewhere in (0, Inf).
			if (_grad_nu_min > 0 && _grad_nu_max < 0) {
				// In this case, the function reaches its mode within [nu_min, nu_max].
				_x_min = _nu_min;
				_x_max = DirectSampling::find_root(_nu_min, _nu_max, grad, _tol);
				_shape = mode;
			} else if (_grad_nu_min > 0 && _grad_nu_max > 0) {
				// In this case, the function is increasing on [nu_min, nu_max].
				_x_min = _nu_min;
				_x_max = _nu_max;
				_shape = increasing;
			} else if (_grad_nu_min < 0 && _grad_nu_max < 0) {
				// In this case, the function is decreasing on [nu_min, nu_max].
				_x_min = _nu_max;
				_x_max = _nu_min;
				_shape = decreasing;
			} else {
				Rcpp::stop("We should not have gotten here");
			}
		} else {
			// In this case, log w(x) is strictly increasing on (0, Inf).
			_x_min = _nu_min;
			_x_max = _nu_max;
			_shape = increasing;
		}

		_log_min = eval(_x_min, true);
		_log_max = eval(_x_max, true);
	}

	// Evaluate the function
	double eval(double x, bool log = false) const
	{
		double out = R_NegInf;
		if (x >= _nu_min && x <= _nu_max) {
			out = _n * x/2 * std::log(x/2) - _n * lgamma(x/2) - _A * x;
		}
		return log ? out :exp(out);
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const
	{
		double x1;
		double x2;
		if (log_a < _log_min || std::isinf(log_a)) {
			return std::pair<double,double>(_nu_min, _nu_max);
		}

		// Given log_a is larger than log_max, so return a nonsensical (empty)
		// interval.
		if (log_a > _log_max) {
			return std::pair<double,double>(_nu_max, _nu_min);
		}

		const DirectSampling::dfd_t& f = [&](double x) -> double {
			double out = R_NegInf;
			if (x >= _nu_min && x <= _nu_max) {
				out = _n * x/2 * log(x/2) - _n * lgamma(x/2) - _A * x - log_a;
			}
			return out;
		};

		if (_shape == mode) {
			// Here there will be two roots which aren't necessarily nu_min or nu_max
			if (eval(_nu_min, true) < log_a) {
				x1 = DirectSampling::find_root(_nu_min, _x_max, f, _tol);
			} else {
				x1 = _nu_min;
			}

			if (eval(_nu_max, true) < log_a) {
				x2 = DirectSampling::find_root(_x_max, _nu_max, f, _tol);
			} else {
				x2 = _nu_max;
			}
		} else if (_shape == increasing) {
			x1 = DirectSampling::find_root(_nu_min, _nu_max, f, _tol);
			x2 = _nu_max;
		} else if (_shape == decreasing) {
			x1 = _nu_min;
			x2 = DirectSampling::find_root(_nu_min, _nu_max, f, _tol);
		} else {
			Rcpp::stop("Unrecognized value of shape");
		}

		// Rprintf("TDistDFWeightFunction::roots(%g) = (%g,%g)\n", log_a, x1, x2);
		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const
	{
		return _log_max;
	}

private:
	double _n;
	double _A;
	double _nu_min;
	double _nu_max;
	double _tol;
	double _grad_nu_min;
	double _grad_nu_max;
	double _x_min;
	double _x_max;
	double _log_min;
	double _log_max;
	TDistDFShape _shape;
};

#endif
