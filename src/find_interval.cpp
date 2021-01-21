#include "find_interval.h"
#include "functionals.h"
#include "bisection.h"

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
//' @export
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
