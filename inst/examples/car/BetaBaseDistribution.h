// [[Rcpp::depends(DirectSampling)]]
#ifndef BETA_BASE_DISTRIBUTION
#define BETA_BASE_DISTRIBUTION

#include "DirectSampling.h"

class BetaBaseDistribution : public DirectSampling::BaseDistribution
{
public:
	BetaBaseDistribution(double a, double b)
	: DirectSampling::BaseDistribution(), _a(a), _b(b)
	{
	}

	// Evaluate the Beta(a, b) density function
	double density(double x, bool log = false) const
	{
		return R::dbeta(x, _a, _b, log);
	}

	// Compute Pr(x1 < X < x2) probability where X ~ Beta(a, b)
	double pr_interval(double x1, double x2, bool log) const
	{
		double lp1 = R::pbeta(x1, _a, _b, true, true);
		double lp2 = R::pbeta(x2, _a, _b, true, true);
		double out = DirectSampling::log_sub2_exp(lp2, lp1);
		return log ? out : std::exp(out);
	}

	// Quantile function of Beta(a, b) truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const
	{
		double p_min = R::pbeta(x_min, _a, _b, false, false);
		double p_max = R::pbeta(x_max, _a, _b, false, false);
		double x = R::qbeta((p_max - p_min)*p + p_min, _a, _b, false, false);
		return std::max(x_min, std::min(x, x_max));
	}

private:
	double _a;
	double _b;
};

#endif

