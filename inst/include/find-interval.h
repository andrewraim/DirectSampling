#ifndef FIND_INTERVAL_H
#define FIND_INTERVAL_H

#include <Rcpp.h>
#include "functionals.h"
#include "bisection.h"

namespace DirectSampling {

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

//' Find Interval
//' 
//' @param x A number 
//' @param cutpoints A sorted vector
//' 
//' @details
//' Elements of cutpoints represent endpoints
//' of adjacent intervals \eqn{[c_0,c_1)}, \eqn{[c_1,c_2)} ...,
//' \eqn{[c_{k},c_{k+1})}. Return the index \code{i} such that \code{x} is in 
//' \eqn{[c_i,c_{i+1})}; or return \code{-1} if \eqn{x < c_0}
//' or \code{k+1} if \eqn{x > c_{k+1}}
//'
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

}

#endif
