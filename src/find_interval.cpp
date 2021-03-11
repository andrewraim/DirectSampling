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

