#include "dgeom.h"

/*
* Here we make the argument is an integer. This is not ideal for general use,
* because, for example, d_dgeom(3.5, 0.5, false) == d_dgeom(3, 0.5, false).
* However, testing whether double exactly represents an integer does not
* seem trivial. See R's built-in dgeom function, for example.
*/
double d_dgeom(int x, double p, bool log)
{
	double out = std::log(p) - std::log(2-p) + abs(x) * std::log(1-p);
	return log ? out : std::exp(out);
}

double p_dgeom(double x, double p)
{
	double xx = floor(x);
	double out;
	if (xx >= 0) {
		out = (1-p) / (2-p) + (1 - pow(1-p, xx+1)) / (2-p);
	} else {
		out =  pow(1-p, -xx) / (2-p);
	}
	return out;
}

int r_dgeom(double p)
{
	return R::rgeom(p) - R::rgeom(p);
}

double q_dgeom(double q, double p)
{
	double x_neg = ceil( -log((2-p)*q) / log(1-p) );
	double x_pos = ceil( (log(1-q) + log(2-p)) / log(1-p) - 1 );
	return  x_neg*(x_neg < 0) + x_pos*(x_pos >= 0);
}

Rcpp::NumericVector d_dgeom(const Rcpp::NumericVector& x, double p, bool log = false)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_dgeom(x(i), p, log);
	}

	return out;
}

Rcpp::NumericVector p_dgeom(const Rcpp::NumericVector& x, double p)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_dgeom(x(i), p);
	}

	return out;
}

Rcpp::NumericVector r_dgeom(unsigned int n, double p)
{
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_dgeom(p);
	}
	return out;
}

Rcpp::NumericVector q_dgeom(const Rcpp::NumericVector& q, double p)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_dgeom(q(i), p);
	}

	return out;
}
