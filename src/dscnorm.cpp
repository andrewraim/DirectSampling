#include <Rcpp.h>
#include "dgeom.h"

double d_dscnorm_unnormalized(int x, double sigma)
{
	double sigma2 = pow(sigma, 2.0);
	int x2 = x * x;
	return exp( -x2/(2*sigma2) );
}

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

Rcpp::NumericVector d_dscnorm_unnormalized(const Rcpp::IntegerVector& x, double sigma)
{
	unsigned n = x.size();
	Rcpp::NumericVector out(n);
	
	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_dscnorm_unnormalized(x(i), sigma);
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
