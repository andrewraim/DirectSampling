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
	double density(double x, bool take_log = false) const {
		return R::dunif(x, _a, _b, take_log);
	}
	// Compute Pr(x1 < X < x2) probability where X ~ Uniform(a, b)
	double pr_interval(double x1, double x2, bool take_log) const {
		double out = R::punif(x2, _a, _b, true, false) - R::punif(x1, _a, _b, true, false);
		if (take_log) { return log(out); } else { return out; }
	}
	// Quantile function of Uniform truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = R::punif(x_min, _a, _b, true, false);
		double p_max = R::punif(x_max, _a, _b, true, false);
		double p_trunc = (p_max - p_min)*p + p_min;
		double x = R::qunif(p_trunc, _a, _b, true, false);
		//Rprintf("p = %g, x_min = %g, x_max = %g, p_min = %g, p_max = %g, p_trunc = %g, x = %g\n",
		//	p, x_min, x_max, p_min, p_max, p_trunc, x);
		return std::max(x_min, std::min(x, x_max));
	}
private:
	double _a;
	double _b;
};

#endif

