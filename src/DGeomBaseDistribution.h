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
	// Compute Pr(x1 < X <= x2) probability where X ~ DGeom(rho)
	double pr_interval(double x1, double x2) const {
		return p_dgeom(x2, _rho) - p_dgeom(x1, _rho);
	}
	// Quantile function of DGeom truncated to [x_min, x_max]
	double q_truncated(double p, double x_min, double x_max) const {
		double p_min = p_dgeom(ceil(x_min) - 1, _rho);
		double p_max = p_dgeom(floor(x_max), _rho);
		double x = q_dgeom((p_max - p_min)*p + p_min, _rho);
		return std::max(ceil(x_min), std::min(x, floor(x_max)));
	}
private:
	double _rho;
};

#endif
