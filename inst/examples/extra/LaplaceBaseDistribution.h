#ifndef LAPLACE_BASE_DISTRIBUTION
#define LAPLACE_BASE_DISTRIBUTION

#include "BaseDistribution.h"
#include "laplace.h"

class LaplaceBaseDistribution : public BaseDistribution
{
public:
	LaplaceBaseDistribution(double lambda)
		: BaseDistribution(), _lambda(lambda)
	{
	}

	// Evaluate the Laplace(0, lambda) density function
	double density(double x, bool log = false) const {
		return d_laplace(x, 0, _lambda, log);
	}

	// Compute Pr(x1 < X < x2) probability where X ~ Laplace(0, lambda)
	double pr_interval(double x1, double x2, double log) const {
		double lp1 = p_laplace(x1, 0, _lambda);
		double lp2 = p_laplace(x2, 0, _lambda);
		double out = DirectSampling::log_sub2_exp(lp2, lp1);
		return log ? out : std::exp(out);
	}

	// Quantile function of Laplace(0, lambda) truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = p_laplace(x_min, 0, _lambda);
		double p_max = p_laplace(x_max, 0, _lambda);
		double x = q_laplace((p_max - p_min)*p + p_min, 0, _lambda);
		return std::max(x_min, std::min(x, x_max));
	}

private:
	double _lambda;
};

#endif

