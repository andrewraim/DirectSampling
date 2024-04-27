// [[Rcpp::depends(DirectSampling)]]
#include "DirectSampling.h"
#include "UniformBaseDistribution.h"
#include "TDistDFWeightFunction.h"

// [[Rcpp::export]]
Rcpp::List r_tdist_unif(unsigned int n, 
	unsigned int m, double A, double nu_min, double nu_max,
	unsigned int N, unsigned int max_rejections, double tol = 1e-10)
{
	TDistDFWeightFunction w(m, A, nu_min, nu_max, 1e-10);
	UniformBaseDistribution g(nu_min, nu_max);
	return DirectSampling::direct_sampler_ar(n, w, g, tol, N, max_rejections,
		"small_rects", 0.5);
}

