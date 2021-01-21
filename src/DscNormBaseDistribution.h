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
	// Compute Pr(x1 < X <= x2) probability where X ~ DscNorm(0, sigma2)
	double pr_interval(double x1, double x2) const {
		return p_dscnorm(x2, _sigma, _tol) - p_dscnorm(x1, _sigma, _tol);
	}
	// Quantile function of DscNorm truncated to [x_min, x_max]
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = p_dscnorm(ceil(x_min) - 1, _sigma, _tol);
		double p_max = p_dscnorm(floor(x_max), _sigma, _tol);
		double x = q_dscnorm((p_max - p_min)*p + p_min, _sigma, _tol);
		return std::max(ceil(x_min), std::min(x, floor(x_max)));
	}
private:
	double _sigma;
	double _tol;
};

#endif
