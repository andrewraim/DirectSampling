#ifndef CAR_WEIGHT_FUNCTION
#define CAR_WEIGHT_FUNCTION

#include <Rcpp.h>
#include "../Rcpp/WeightFunction.h"
#include "../Rcpp/find_root.h"

/*
* This class is used in root-finding
*/
class CAREval : public Functional1
{
public:
	CAREval(const Rcpp::NumericVector& lambda, double C, double x_lower, double x_upper, double log_a)
		: _lambda(lambda), _C(C), _x_lower(x_lower), _x_upper(x_upper), _log_a(log_a)
	{
	}

	double operator()(double x) const {
		double out = -INFINITY;
		if (x >= _x_lower && x <= _x_upper) {
			out = _C*x + 0.5 * Rcpp::sum(Rcpp::log((1 - _lambda*x))) - _log_a;
		}
		return out;
	}

private:
	Rcpp::NumericVector _lambda;
	double _C;
	double _x_lower;
	double _x_upper;
	double _log_a;
};

/*
* This class is used in root-finding
*/
class CARGrad : public Functional1
{
public:
	CARGrad(const Rcpp::NumericVector& lambda, double C, double x_lower, double x_upper)
		: _lambda(lambda), _C(C), _x_lower(x_lower), _x_upper(x_upper)
	{
	}

	double operator()(double x) const {
		double out = 0;
		if (x >= _x_lower && x <= _x_upper) {
			out = _C - 0.5 * Rcpp::sum(_lambda / (1 - _lambda*x));
		}
		return out;
	}

private:
	Rcpp::NumericVector _lambda;
	double _C;
	double _x_lower;
	double _x_upper;
};

class CARWeightFunction : public WeightFunction
{
public:
	CARWeightFunction(const Rcpp::NumericVector& lambda, double C,
		double x_lower, double x_upper, double tol)
		: WeightFunction(), _lambda(lambda), _C(C), _x_lower(x_lower),
			_x_upper(x_upper), _x_min(x_lower), _x_max(x_upper), _log_min(-INFINITY),
			_log_max(INFINITY), _tol(tol)
	{
		init();
	}

	void init() {
		CARGrad grad(_lambda, _C, _x_lower, _x_upper);

		// x_min occurs at one of the two endpoints.
		// x_max occurs between the two endpoints if grad(0) > 0 and grad(1) < 0.
		// Otherwise x_max is one of the endpoints.
		double grad_lower = grad(_x_lower);
		double grad_upper = grad(_x_upper);
		if (grad_lower > 0 && grad_upper < 0) {
			// Recall that grad is a decreasing function.
			// In this case, there is a max strictly between (x_lower, x_upper)
			// The min occurs at one of the endpoints
			_x_max = find_root(_x_lower, _x_upper, grad, _tol);
			_log_max = eval(_x_max, true);

			double log_lower = eval(_x_lower, true);
			double log_upper = eval(_x_upper, true);
			if (log_lower < log_upper) {
				_x_min = _x_lower;
				_log_min = log_lower;
			} else {
				_x_min = _x_upper;
				_log_min = log_upper;
			}
		} else if (grad_lower < 0) {
			// In this case, the gradient is negative throughout [x_lower, x_upper]
			// so log w(x) is decreasing.
			_x_min = _x_upper;
			_x_max = _x_lower;
			_log_max = eval(_x_max, true);
			_log_min = eval(_x_min, true);
		} else if (grad_upper > 0) {
			// In this case, the gradient is positive throughout [x_lower, x_upper]
			// so log w(x) is increasing.
			_x_min = _x_lower;
			_x_max = _x_upper;
			_log_max = eval(_x_max, true);
			_log_min = eval(_x_min, true);
		} else {
			Rcpp::stop("I didnt program this case yet");
		}
	}

	// Evaluate the function
	double eval(double x, bool take_log = false) const {
		double out = -INFINITY;
		if (x >= _x_lower && x <= _x_upper) {
			out = _C*x + 0.5 * Rcpp::sum(Rcpp::log((1 - _lambda*x)));
		}
		if (take_log) { return out; } else { return exp(out); }
	}

	// Return the roots of the equation w(x) = a, which is equivalent to
	// log(w(x)) = log(a). Roots are returned in increasing order.
	std::pair<double,double> roots(double log_a) const {
		double x1;
		double x2;
		if (log_a < _log_min || std::isinf(log_a)) {
			return std::pair<double,double>(_x_lower, _x_upper);
		}

		if (log_a < _log_min) {
			return std::pair<double,double>(_x_lower, _x_upper);
		}

		// Given log_a is larger than log_max, so return a nonsensical (empty)
		//interval.
		if (log_a > _log_max) {
			return std::pair<double,double>(_x_upper, _x_lower);
		}

		CAREval f(_lambda, _C, _x_lower, _x_upper, log_a);

		if (f(_x_lower) < 0) {
			x1 = find_root(_x_lower, _x_max, f, _tol);
		} else {
			x1 = _x_lower;
		}

		if (f(_x_upper) < 0) {
			x2 = find_root(_x_max, _x_upper, f, _tol);
		} else {
			x2 = _x_upper;
		}

		return std::pair<double,double>(x1, x2);
	}

	// log_c is the maximum value of the function log[w(x)]
	double log_c() const {
		return _log_max;
	}
private:
	Rcpp::NumericVector _lambda;
	double _C;
	double _x_lower;
	double _x_upper;
	double _x_min;
	double _x_max;
	double _log_min;
	double _log_max;
	double _tol;
};

#endif

