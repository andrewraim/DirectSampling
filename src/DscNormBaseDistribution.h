#ifndef DSCNORM_BASE_DISTRIBUTION
#define DSCNORM_BASE_DISTRIBUTION

#include "BaseDistribution.h"
#include "dscnorm.h"

class DscNormBaseDistribution : public BaseDistribution
{
public:
	DscNormBaseDistribution(double sigma)
		: BaseDistribution(), _sigma(sigma), _tol(1e-10)
	{
	}
	DscNormBaseDistribution(double sigma, double tol)
		: BaseDistribution(), _sigma(sigma), _tol(tol)
	{
	}
	// Evaluate the DscNorm(0, sigma2) density function
	double density(double x, bool take_log = false) const {
		return d_dscnorm(x, _sigma, _tol, take_log, true);
	}
	// Compute Pr(x1 < X < x2) probability where X ~ DscNorm(0, sigma2)
	// This calculation should exclude endpoints if they are integers
	double pr_interval(double x1, double x2) const {
		double a = floor(x1 + 1);
		double b = ceil(x2 - 1);
		return p_dscnorm(b, _sigma, _tol) - p_dscnorm(a - 1, _sigma, _tol);
	}
	// Quantile function of DscNorm truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		double a = floor(x_min + 1);
		double b = ceil(x_max - 1);
		double p_min = p_dscnorm(a - 1, _sigma, _tol);
		double p_max = p_dscnorm(b, _sigma, _tol);
		double x = q_dscnorm((p_max - p_min)*p + p_min, _sigma, _tol);
		return std::max(ceil(x_min), std::min(x, floor(x_max)));
	}
private:
	double _sigma;
	double _tol;
};

#endif
