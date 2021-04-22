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

