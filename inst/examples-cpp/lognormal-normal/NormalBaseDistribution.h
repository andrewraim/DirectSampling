// [[Rcpp::depends(DirectSampling)]]
#ifndef NORMAL_BASE_DISTRIBUTION
#define NORMAL_BASE_DISTRIBUTION

#include "DirectSampling.h"

class NormalBaseDistribution : public DirectSampling::BaseDistribution
{
public:
	NormalBaseDistribution(double sigma)
		: DirectSampling::BaseDistribution(), _sigma(sigma)
	{
	}
	// Evaluate the N(0, sigma2) density function
	double density(double x, bool take_log = false) const {
		return R::dnorm(x, 0, _sigma, take_log);
	}
	// Compute Pr(x1 < X < x2) probability where X ~ N(0, sigma2)
	double pr_interval(double x1, double x2, bool take_log) const {
		double out = R::pnorm(x2, 0, _sigma, true, false) - R::pnorm(x1, 0, _sigma, true, false);
		return take_log ? log(out) : out;
	}
	// Quantile function of N truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = R::pnorm(x_min, 0, _sigma, true, false);
		double p_max = R::pnorm(x_max, 0, _sigma, true, false);
		double p_trunc = (p_max - p_min)*p + p_min;
		double x = R::qnorm(p_trunc, 0, _sigma, true, false);
		return std::max(x_min, std::min(x, x_max));
	}
private:
	double _sigma;
};

#endif
