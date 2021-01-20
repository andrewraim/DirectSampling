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

// Currently assumes both roots are real
// [[Rcpp::export]]
Rcpp::NumericVector quadratic_roots(double A, double B, double C)
{
	return Rcpp::NumericVector::create(
		(-B + sqrt(B*B - 4*A*C)) / (2*A),
		(-B - sqrt(B*B - 4*A*C)) / (2*A)
	);
}

Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease)
{
	Rcpp::NumericVector sorted = clone(x).sort(decrease);
	return Rcpp::match(sorted, x);
}

// Compute log(p[0] + ... + p[k-1]) given log(p[0]), ..., log(p[k-1]).
// Treat p[0] as the "normalizing" probability.
// Uses the method from
// https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation/subtraction
// for stable calculation on the log scale.
// [[Rcpp::export]]
double logsumprobs(const Rcpp::NumericVector& logprob)
{
	unsigned int k = logprob.size();
	const Rcpp::IntegerVector& idx = Rcpp::seq_len(k-1);
	const Rcpp::NumericVector& a = logprob[idx];
	const Rcpp::NumericVector& b = Rcpp::rep(logprob(0), idx.length());
	double s = Rcpp::sum(Rcpp::exp(a - b));
	return logprob(0) + log1p(s);
}
