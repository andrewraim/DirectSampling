#ifndef NORMAL_BASE_DISTRIBUTION
#define NORMAL_BASE_DISTRIBUTION

#include "BaseDistribution.h"

class NormalBaseDistribution : public BaseDistribution
{
public:
	NormalBaseDistribution(double sigma)
		: BaseDistribution(), _sigma(sigma)
	{
	}
	// Evaluate the N(0, sigma2) density function
	double density(double x, bool take_log = false) const {
		return R::dnorm(x, 0, _sigma, take_log);
	}
	// Compute Pr(x1 < X <= x2) probability where X ~ N(0, sigma2)
	double pr_interval(double x1, double x2) const {
		return R::pnorm(x2, 0, _sigma, true, false) - R::pnorm(x1, 0, _sigma, true, false);
	}
	// Quantile function of N truncated to [x_min, x_max]
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = R::pnorm(x_min, 0, _sigma, true, false);
		double p_max = R::pnorm(x_max, 0, _sigma, true, false);
		double x = R::qnorm((p_max - p_min)*p + p_min, 0, _sigma, true, false);
		return std::max(x_min, std::min(x, x_max));
	}
private:
	double _sigma;
};

#endif

