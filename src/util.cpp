#include "util.h"

Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
	Rcpp::NumericVector z(x.size() + y.size());
	for (unsigned int i = 0; i < x.size(); i++) {
		z(i) = x(i);
	}
	for (unsigned int i = 0; i < y.size(); i++) {
		z(i + x.size()) = y(i);
	}
	return z;
}

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x)
{
	std::vector<unsigned int> idx;
	for (unsigned int i = 0; i < x.size(); i++) {
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
	Rcpp::NumericVector out(x.size());
	for (unsigned int i = 0; i < out.size(); i++) {
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
