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

double bisection(double x_lo, double x_hi,
	const Predicate& pred,
	const Functional2& mid,
	const Functional2& dist, double tol)
{
	double x = mid(x_lo, x_hi);

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		bool ind = pred(x);
		x_lo = ind*x_lo + (1-ind)*x;
		x_hi = ind*x + (1-ind)*x_hi;
		x = mid(x_lo, x_hi);
	}

	if (dist(x_lo, x_hi) > tol && (x <= x_lo || x >= x_hi)) {
		Rcpp::stop("Numerical overflow in bisection. tol may be too small");
	}

	return x;
}

class FloorMidpoint : public Functional2
{
public:
	virtual double operator()(double x, double y) const
	{
		return floor((y + x) / 2);
	}
};

class FindIntervalPredicate : public Predicate
{
public:
	FindIntervalPredicate(const Rcpp::NumericVector& cutpoints, double x)
		: _cutpoints(cutpoints), _x(x)
	{
	}
	bool operator()(double idx) const {
		return _cutpoints(int(idx)) >= _x;
	}
private:
	const Rcpp::NumericVector& _cutpoints;
	double _x;
};

// Assume cutpoints is a sorted vector, where elements represent endpoints
// of adjacent intervals [c_0,c_1), [c_1,c_2) ..., [c_{k},c_{k+1}). Return
// the index i such that x is in [c_i,c_{i+1}); or return -1 if x < c_0
// or k+1 if x > c_{k+1}
// [[Rcpp::export]]
int find_interval(double x, const Rcpp::NumericVector& cutpoints)
{
	unsigned int N = cutpoints.size();
	unsigned int k = N-1;

	if (x < cutpoints(0)) {
		return -1;
	} else if (x > cutpoints(k)) {
		return N;
	}

	FindIntervalPredicate pred(cutpoints, x);
	FloorMidpoint mid;
	IntervalLength dist;
	double out = bisection(0, k, pred, mid, dist, 1);
	return int(out);
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

// Try to make this version as fast as possible. E.g. avoid
// calling the one above.
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_mat(const Rcpp::NumericMatrix& P)
{
	unsigned int n = P.rows();
	unsigned int k = P.cols();

	const Rcpp::NumericVector& u = Rcpp::runif(n);
	Rcpp::IntegerVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		const Rcpp::NumericVector& cp = Rcpp::cumsum(P.row(i));
		if (cp(k-1) <= 0) {
			// Watch out for the case where all probs are zero
			z(i) = NA_INTEGER;
		} else {
			z(i) = (u(i) < 1)*Rcpp::sum(u(i) > cp) + (u(i) >= 1)*(k-1);
		}
	}

	return z;
}

// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_logp(unsigned int n, const Rcpp::NumericVector& logprob)
{
	const Rcpp::IntegerVector& idx_sort = order(logprob, true);
	const Rcpp::NumericVector& lp = logprob[idx_sort];

	const Rcpp::NumericVector& logu = Rcpp::log(Rcpp::runif(n));
	Rcpp::IntegerVector idx_draw(n);

	unsigned int k = lp.size();
	Rcpp::NumericVector log_cp(k);
	log_cp(0) = lp(0);
	for (unsigned int j = 1; j < k; j++) {
		Rcpp::IntegerVector idx = Rcpp::seq(0, j);
		log_cp(j) = logsumprobs(lp[idx]);
	}

	for (unsigned int i = 0; i < n; i++) {
		idx_draw(i) = (logu(i) < 0)*Rcpp::sum(logu(i) > log_cp) + (logu(i) >= 0)*(k-1);
	}
	return idx_sort[idx_draw];
}

Rcpp::IntegerVector rcateg_logp_mat(const Rcpp::NumericMatrix& logprob)
{
	unsigned int n = logprob.rows();
	Rcpp::IntegerVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::IntegerVector res = rcateg_logp(1, logprob.row(i));
		z(i) = res(0);
	}

	return z;
}



unsigned int rcateg(const Rcpp::NumericVector& p)
{
	unsigned int k = p.size();
	double u = R::runif(0, 1);
	const Rcpp::NumericVector& cp = Rcpp::cumsum(p);

	if (cp(k-1) <= 0) {
		// Watch out for the case where all probs are zero
		return NA_INTEGER;
	}

	// This is a sneaky way of finding the max index such that u(i)> cp.
	// It is possible for the RNG to draw exactly 1; we have to handle this
	// case specially.
	return (u < 1)*Rcpp::sum(u > cp) + (u >= 1)*(k-1);
}

// [[Rcpp::export]]
Rcpp::IntegerVector rcateg(unsigned int n, const Rcpp::NumericVector& p)
{
	unsigned int k = p.size();
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	const Rcpp::NumericVector& cp = Rcpp::cumsum(p);

	if (cp(k-1) <= 0) {
		// Watch out for the case where all probs are zero
		return Rcpp::rep(NA_INTEGER, n);
	}
	
	Rcpp::IntegerVector idx(n);
	for (unsigned int i = 0; i < n; i++) {
		// This is a sneaky way of finding the max index such that u(i)> cp.
		// It is possible for the RNG to draw exactly 1; we have to handle this
		// case specially.
		idx(i) = (u(i) < 1)*Rcpp::sum(u(i) > cp) + (u(i) >= 1)*(k-1);
	}
	return idx;
}
