#include "util.h"

Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
	unsigned int m = x.size();
	unsigned int n = y.size();

	Rcpp::NumericVector z(m + n);
	for (unsigned int i = 0; i < m; i++) {
		z(i) = x(i);
	}
	for (unsigned int i = 0; i < n; i++) {
		z(i + m) = y(i);
	}
	return z;
}

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x)
{
	unsigned int n = x.size();
	std::vector<unsigned int> idx;

	for (unsigned int i = 0; i < n; i++) {
		if (x(i)) {
			idx.push_back(i);
		}
	}

	Rcpp::IntegerVector out(idx.begin(), idx.end());
	return out;
}

// For some reason, Rcpp/Sugar cumsum does not compile for me...
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x)
{
	double sum = 0;
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		sum += x(i);
		out(i) = sum;
	}
	return out;
}

Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease)
{
	Rcpp::NumericVector sorted = clone(x).sort(decrease);
	return Rcpp::match(sorted, x);
}

Rcpp::IntegerVector q_discrete(const Rcpp::NumericVector& q, const Rcpp::NumericVector& cp)
{
	unsigned int n = q.size();
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_discrete(q(i), cp);
	}

	return out;
}

unsigned int q_discrete(double q, const Rcpp::NumericVector& cp)
{
	unsigned int k = cp.size();
	if (q > cp(k-1)) {
		Rcpp::stop("q > max(cp)");
	}

	// Otherwise do a binary search
	unsigned int x_lo = 0;
	unsigned int x_hi = k-1;
	unsigned int x = (unsigned int)(floor((x_hi + x_lo) / 2.0));
	while (x_hi - x_lo > 1) {
		bool ind = (cp(x) >= q);
		x_lo = (1 - ind) * x + ind * x_lo;
		x_hi = ind * x + (1 - ind) * x_hi;
		x = (unsigned int)(floor((x_hi + x_lo) / 2.0));
	}

	if (cp(x_lo) >= q) {
		return x_lo;
	} else {
		return x_hi;
	}
}
