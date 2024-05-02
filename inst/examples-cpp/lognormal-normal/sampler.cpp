// [[Rcpp::depends(DirectSampling)]]
#include "DirectSampling.h"
#include "LognormalWeightFunction.h"
#include "NormalBaseDistribution.h"

//' Lognormal-Normal Direct Sampler
//' 
//' Lognormal-Normal direct sampler with accept-reject implemented in C++,
//' 
//' @param n Number of draws to generate.
//' @param z A vector of points to evaluate.
//' @param mu Mean parameter for weight function.
//' @param sigma2 standard deviation parameter for weight function.
//' @param lambda Standard deviation parameter for base distribution.
//' @param N Number of knots to use in approximation for \eqn{p(u)}.
//' @param max_rejections Maximum number of allowed rejections in accept-reject.
//' 
//' @details
//' Implementation of the a direct sampler with weight function
//' \deqn{
//' w(\xi) = \frac{1}{z - \xi} \exp\left\{ -\frac{1}{2\sigma^2}
//' (\log (z - \xi) - \mu)^2 \right\} \mathrm{I}\{\xi < z \}
//' }
//' corresponding to Lognormal density and base distribution
//' \deqn{
//' g(\xi) = N(\xi \mid 0, \lambda^2).
//' }
//' 
//' @returns
//' A list with element \code{x} a vector of the accepted draws and
//' element \code{rejections} a vector of counts recording the number of
//' rejections required to obtain corresponding draws in \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu,
	double sigma2, double lambda, unsigned int N, unsigned int max_rejections)
{
	LognormalWeightFunction w(z, mu, sigma2);
	NormalBaseDistribution g(lambda);

	// These typically won't need to be changed
	double tol = 1e-8;
	std::string method = "small_rects";
	double priority_weight = 0.5;

	return DirectSampling::direct_sampler_ar(n, w, g, tol, N,
		max_rejections, method, priority_weight);
}

