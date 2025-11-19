#ifndef DSCNORM_BASE_DISTRIBUTION
#define DSCNORM_BASE_DISTRIBUTION

#include "BaseDistribution.h"
#include "dscnorm.h"
#include "util.h"

class DscNormBaseDistribution : public BaseDistribution
{
public:
	DscNormBaseDistribution(double sigma)
	: BaseDistribution(), _sigma(sigma), _tol(1e-10)
	{
		init();
	}

	DscNormBaseDistribution(double sigma, double tol)
	: BaseDistribution(), _sigma(sigma), _tol(tol)
	{
		init();
	}

	// Calling p_dscnorm and q_dscnorm directly results in cumulative
	// probabilities being recomputed many times. We precompute those quantities
	// here and reuse them to try to avoid this extra work.
	void init() {
		// Set min and max values based on truncated support
		_x_hi = hi_dscnorm(_sigma, _tol);
		_x_lo = -_x_hi;
		unsigned int n = _x_hi - _x_lo + 1;

		_cumprobs = Rcpp::NumericVector(n);

		// First set _cumprobs to unnormalized probabilities.
		double cp_unnorm = 0;
		double sigma2 = _sigma * _sigma;
		for (unsigned int i = 0; i < n; i++) {
			int v = _x_lo + i;
			double fx = exp( -v*v / (2*sigma2) );
			cp_unnorm += fx;
			_cumprobs(i) = cp_unnorm;
		}

		// Now divide _cumprobs by the normalizing constant.
		_cumprobs = _cumprobs / cp_unnorm;
	}

	// Evaluate the DscNorm(0, sigma2) density function
	double density(double x, bool log = false) const {
		return d_dscnorm(x, _sigma, _tol, log, true);
	}

	// Compute Pr(x1 < X < x2) probability where X ~ DscNorm(0, sigma2)
	// This calculation should exclude endpoints if they are integers
	double pr_interval(double x1, double x2, bool log) const {
		double a = floor(x1 + 1);
		double b = ceil(x2 - 1);
		double out =  p_local(b) - p_local(a - 1);
		return log ? out : std::exp(out);
	}

	// Quantile function of DscNorm truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		if (ceil(x_min) > floor(x_max)) {
			Rcpp::stop("There are no integers between x_min = %g and x_max = %g\n", x_min, x_max);
		}
		double a = floor(x_min + 1);
		double b = ceil(x_max - 1);
		double p_min = p_local(a - 1);
		double p_max = p_local(b);
		int x = q_local((p_max - p_min)*p + p_min);
		return std::max(ceil(x_min), std::min(double(x), floor(x_max)));
	}

	// CDF of DscNorm which avoids recomputing cumprobs
	double p_local(double x) const {
		int z = floor(x);
		if (z < _x_lo) {
			return 0;
		} else if (z > _x_hi) {
			return 1;
		} else {
			unsigned int idx = z - _x_lo;
			return _cumprobs(idx);
		}
	}

	// Quantile function of DscNorm which avoids recomputing cumprobs
	int q_local(double p) const {
		unsigned int idx = q_discrete(p, _cumprobs);
		return idx + _x_lo;
	}

private:
	double _sigma;
	double _tol;
	int _x_lo;
	int _x_hi;
	Rcpp::NumericVector _cumprobs;
};

#endif
