// [[Rcpp::depends(DirectSampling, fntl)]]
#include "DirectSampling.h"
#include "UniformBaseDistribution.h"
#include "TDistDFWeightFunction.h"

namespace ds = DirectSampling;

// [[Rcpp::export]]
Rcpp::List r_tdist_unif(unsigned int n, 
	unsigned int m, double A, double nu_min, double nu_max,
	unsigned int N, unsigned int max_rejects, double tol = 1e-10)
{
	TDistDFWeightFunction w(m, A, nu_min, nu_max, 1e-10);
	UniformBaseDistribution g(nu_min, nu_max);

	ds::stepdown_args args1;
	args1.N = N;
	args1.tol = tol;
	ds::Stepdown step(w, g, args1);

	ds::rejection_args args2;
	args2.max_rejects = max_rejects;

	const Rcpp::List& out = ds::rejection(step, n, args2);
	return out;
}

