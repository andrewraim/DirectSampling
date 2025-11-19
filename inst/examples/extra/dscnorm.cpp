#include "dscnorm.h"
#include "dgeom.h"
#include "find_interval.h"

int hi_dscnorm(double sigma, double tol)
{
	return ceil(R::qnorm(tol / 2, 0, sigma, false, false));
}

double d_dscnorm(int x, double sigma, double tol, bool log, bool normalize)
{
	double sigma2 = sigma * sigma;
	int x_hi = hi_dscnorm(sigma, tol);
	int x_lo = -x_hi;

	double norm_const = 1;
	if (normalize) {
		norm_const = 0;
		for (int v = x_lo; v <= x_hi; v++) {
			norm_const += exp( -v*v / (2*sigma2) );
		}
	}

	double out = exp( -x*x / (2*sigma2) ) / norm_const;
	return log ? std::log(out) : out;
}

double p_dscnorm(int x, double sigma, double tol)
{
	double sigma2 = sigma * sigma;
	int x_hi = hi_dscnorm(sigma, tol);
	int x_lo = -x_hi;
	unsigned int n = x_hi - x_lo + 1;

	if (x < x_lo) {
		return 0;
	} else if (x > x_hi) {
		return 1;
	}

	double norm_const = 0;
	double cp_unnorm = 0;
	for (unsigned int i = 0; i < n; i++) {
		int v = x_lo + i;
		double fx = exp( -v*v / (2*sigma2) );
		norm_const += fx;
		cp_unnorm += fx * (v <= x);
	}

	return cp_unnorm / norm_const;
}

/*
* Here we don't need the tolerance. We use the exact rejection sampling
* algorithm given by Cannone et al <https://arxiv.org/abs/2004.00010>
*/
int r_dscnorm(double sigma)
{
	double sigma2 = pow(sigma, 2.0);
	double t = floor(sigma) + 1;
	double rho = 1 - exp(-1/t);
	int y;
	double u;

	bool accept = false;
	while (!accept) {
		y = r_dgeom(rho);
		u = R::runif(0, 1);
		accept = u < exp(-1/(2*sigma2) * pow(abs(y) - sigma2 / t, 2));
	}

	return y;
}

int q_dscnorm(double q, double sigma, double tol)
{
	double sigma2 = sigma * sigma;
	int x_hi = hi_dscnorm(sigma, tol);
	int x_lo = -x_hi;
	unsigned int n = x_hi - x_lo + 1;

	Rcpp::NumericVector probs_unnorm(n);

	for (unsigned int i = 0; i < n; i++) {
		int v = x_lo + i;
		probs_unnorm(i) = exp( -v*v / (2*sigma2) );
	}

	const Rcpp::NumericVector& probs = probs_unnorm / Rcpp::sum(probs_unnorm);

	// Take x to be the first value where the cumulative probability (cp) is at
	// least as large as the target q. This should return x_lo if q = 0 or x_hi
	// if q = 1.
	double cp = 0;
	int x = x_lo;
	for (unsigned int i = 0; i < n && cp < q; i++) {
		cp += probs(i);
		x = x_lo + i;
	}

	return x;
}

/*
* The functions below could be more efficient if they didn't recompute the
* normalizing constant over and over again.
*/

Rcpp::NumericVector d_dscnorm(const Rcpp::NumericVector& x, double sigma,
	double tol, bool log, bool normalize)
{
	unsigned n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_dscnorm(x(i), sigma, tol, log, normalize);
	}

	return out;
}

Rcpp::NumericVector p_dscnorm(const Rcpp::NumericVector& x, double sigma, double tol)
{
	unsigned n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_dscnorm(x(i), sigma, tol);
	}

	return out;
}

Rcpp::IntegerVector r_dscnorm(unsigned int n, double sigma)
{
	Rcpp::IntegerVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_dscnorm(sigma);
	}

	return out;
}

Rcpp::IntegerVector q_dscnorm(const Rcpp::NumericVector& q, double sigma, double tol)
{
	unsigned n = q.size();
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_dscnorm(q(i), sigma, tol);
	}

	return out;
}
