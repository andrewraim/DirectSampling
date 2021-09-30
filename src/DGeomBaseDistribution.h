#ifndef DGEOM_BASE_DISTRIBUTION
#define DGEOM_BASE_DISTRIBUTION

#include "BaseDistribution.h"
#include "dgeom.h"

class DGeomBaseDistribution : public BaseDistribution
{
public:
	DGeomBaseDistribution(double rho)
		: BaseDistribution(), _rho(rho)
	{
	}
	// Evaluate the DGeom(rho) density function
	double density(double x, bool take_log = false) const {
		return d_dgeom(x, _rho, take_log);
	}
	// Compute Pr(x1 < X < x2) probability where X ~ DGeom(rho)
	// This calculation should exclude endpoints if they are integers
	double pr_interval(double x1, double x2) const {
		double a = floor(x1 + 1);
		double b = ceil(x2 - 1);
		return p_dgeom(b, _rho) - p_dgeom(a - 1, _rho);
	}
	// Quantile function of DGeom truncated to (x_min, x_max)
	double q_truncated(double p, double x_min, double x_max) const {
		if (ceil(x_min) > floor(x_max)) {
			Rcpp::stop("There are no integers between x_min = %g and x_max = %g\n", x_min, x_max);
		}
		double a = floor(x_min + 1);
		double b = ceil(x_max - 1);
		double p_min = p_dgeom(a - 1, _rho);
		double p_max = p_dgeom(b, _rho);
		double x = q_dgeom((p_max - p_min)*p + p_min, _rho);
		return std::max(ceil(x_min), std::min(x, floor(x_max)));
	}
private:
	double _rho;
};

#endif
