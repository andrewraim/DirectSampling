// [[Rcpp::depends(DirectSampling)]]
#ifndef UNIFORM_BASE_DISTRIBUTION
#define UNIFORM_BASE_DISTRIBUTION

#include "DirectSampling.h"

class UniformBaseDistribution : public DirectSampling::BaseDistribution
{
public:
	UniformBaseDistribution(double a, double b)
	: DirectSampling::BaseDistribution(), _a(a), _b(b)
	{
	}

	// Evaluate the Uniform(a, b) density function
	double density(double x, bool log = false) const
	{
		return R::dunif(x, _a, _b, log);
	}

	// Compute Pr(x1 < X < x2) probability where X ~ Uniform(a, b)
	double pr_interval(double x1, double x2, bool log) const
	{
		double lp1 = R::punif(x1, _a, _b, true, true);
		double lp2 = R::punif(x2, _a, _b, true, true);
		double out = DirectSampling::log_sub2_exp(lp2, lp1);
		return log ? out : std::exp(out);
	}

	// Quantile function of Uniform truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const
	{
		double p_min = R::punif(x_min, _a, _b, true, false);
		double p_max = R::punif(x_max, _a, _b, true, false);
		double p_trunc = (p_max - p_min)*p + p_min;
		double x = R::qunif(p_trunc, _a, _b, true, false);
		return std::max(x_min, std::min(x, x_max));
	}

private:
	double _a;
	double _b;
};

#endif

