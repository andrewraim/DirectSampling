#ifndef DIRECT_SAMPLING_STEPDOWN_H
#define DIRECT_SAMPLING_STEPDOWN_H

#include <Rcpp.h>
#include <math.h>
#include <queue>
#include "WeightFunction.h"
#include "BaseDistribution.h"
#include "util.h"
#include "log-sum-exp.h"
#include "typedefs.h"
#include "bisection.h"
#include "find-interval.h"
#include "Interval.h"
#include "stepdown-args.h"

namespace DirectSampling {

/*
* This class represents a step function approximation to the non-increasing
* Pr(A_u) for u in [0,1].
*/
class Stepdown {
public:
	Stepdown(const WeightFunction& w, const BaseDistribution& g,
		const stepdown_args& args = stepdown_args());

	Stepdown(const WeightFunction& w, const BaseDistribution& g,
		double tol, unsigned int N, const std::string& method,
		double priority_weight);

	double get_norm_const() const;
	const Rcpp::NumericVector& get_cum_probs() const;
	const Rcpp::NumericVector& get_log_x_vals() const;
	const Rcpp::NumericVector& get_log_h_vals() const;
	const Rcpp::NumericVector& get_knot_order() const;

	const WeightFunction& get_weight() const { return _w; }
	const BaseDistribution& get_base() const { return _g; }

	void add(double log_u);
	double log_p(double log_u) const;
	double density(double log_x, bool log, bool normalize) const;
	double cdf(double x) const;
	double quantile(double p, bool log) const;
	double draw(bool log) const;
	Rcpp::NumericVector draw(unsigned int n, bool log) const;

private:
	void init_equal_steps(double log_L, double log_U, double log_prob_max);
	void init_small_rects(double log_L, double log_U, double log_prob_max);
	void update();

	const WeightFunction& _w;
	const BaseDistribution& _g;
	double _tol;
	unsigned int _N;
	Rcpp::NumericVector _log_x_vals;
	Rcpp::NumericVector _log_h_vals;
	Rcpp::NumericVector _knot_order;
	Rcpp::NumericVector _cum_probs;
	double _norm_const;
	double _priority_weight;
};

inline void print(std::priority_queue<Interval> q)
{
	while(!q.empty()) {
		Interval el = q.top();
		q.pop();
		el.print();
	}
}

inline Stepdown::Stepdown(const WeightFunction& w, const BaseDistribution& g,
	const stepdown_args& args)
: Stepdown(w, g, args.tol, args.N, args.method, args.priority_weight)
{
}

inline Stepdown::Stepdown(const WeightFunction& w, const BaseDistribution& g,
	double tol, unsigned int N, const std::string& method, double priority_weight)
	: _w(w), _g(g), _tol(tol), _N(N), _log_x_vals(), _log_h_vals(), _knot_order(),
	  _cum_probs(), _priority_weight(priority_weight)
{
	double log_prob;
	std::pair<double,double> endpoints;
	
	/* Geometric mean */
	const midpoint_t& midpoint = [](double log_x, double log_y) -> double {
		return 0.5 * (log_x + log_y);
	};

	/* Interval length on the log-scale */
	const distance_t& dist = [](double log_x, double log_y) -> double {
		return log_sub2_exp(log_y, log_x);
	};

	// First, make sure the widest possible A_u intersects with [x_lo, x_hi]. If it
	// doesn't, the rest of the algorithm won't work, so bail out.
	double log_prob_max = log_p(-INFINITY);
	if (std::isinf(log_prob_max)) {
		Rcpp::stop("Could not find any u such that P(X in A_u) > 0");
	}

	// Find L, the largest point where where P(A_L) = prob_max.
	// First find a finite lower bound for log_L.
	double log_L_lo = -1;
	double log_L_hi = 0;
	log_prob = log_p(log_L_lo);
	while (log_prob < log_prob_max) {
		log_L_lo *= 2;
		log_prob = log_p(log_L_lo);
	}

	// Do a bisection search to find log_L between log_L_lo and log_L_hi.
	const predicate_t& predL = [&](double log_u) -> bool {
		return this->log_p(log_u) < log_prob_max;
	};
	double log_delta1 = std::min(log_L_lo, log(tol));
	double log_L = bisection(log_L_lo, log_L_hi, predL, midpoint, dist, log_delta1);

	// Do a bisection search to find U, the smallest point where P(A_U) = 0.
	const predicate_t& predU = [&](double log_u) -> bool {
		return this->log_p(log_u) < log(1e-10) + log_p(log_L);
	};
	double log_delta2 = std::min(log_L, log(tol));
	double log_U = bisection(log_L, 0, predU, midpoint, dist, log_delta2);

	// Now fill in points between L and U
	if (method == "equal_steps") {
		init_equal_steps(log_L, log_U, log_prob_max);
	} else if (method == "small_rects") {
		init_small_rects(log_L, log_U, log_prob_max);
	} else {
		char msg[256];
		sprintf(msg, "Unknown method: %s. Currently support equal_steps and small_rects\n",
			method.c_str());
		Rcpp::stop(msg);
	}

	update();
}

inline void Stepdown::init_equal_steps(double log_L, double log_U, double log_prob_max)
{
	unsigned int N = _N;
	double L = exp(log_L);
	double U = exp(log_U);

	_log_x_vals = Rcpp::NumericVector(N+2);
	_log_h_vals = Rcpp::NumericVector(N+2);

	// Add pieces for the interval [-Inf, 0) where h(u) = 0, and
	// the interval [0, L) where h(u) = 1
	_log_x_vals(0) = -INFINITY;
	_log_h_vals(0) = log_prob_max;

	for (unsigned int i = 0; i < N+1; i++) {
		double prop = i / double(N);
		double log_u = log(L + prop * (U - L));
		_log_x_vals(i+1) = log_u;
		_log_h_vals(i+1) = log_p(log_u);
	}

	_knot_order = Rcpp::seq_len(N+2);
}

inline void Stepdown::init_small_rects(double log_L, double log_U, double log_prob_max)
{
	unsigned int N = _N;
	double pw = _priority_weight;

	/* Geometric mean */
	const midpoint_t& midpoint = [](double log_x, double log_y) -> double {
		return 0.5 * (log_x + log_y);
	};

	// This queue should be in max-heap order by area
	std::priority_queue<Interval> q;
	q.push(Interval(log_L, log_U, log_p(log_L), log_p(log_U), pw));

	/*
	* Try to be efficient by preallocating log_x_vals and
	* log_h_vals to the (maximum) size needed, then copying these temporary
	* structures to _log_x_vals and _log_h_vals afterward. (It looks like
	* Rcpp NumericVectors cannot have their allocation controlled in this way).
	*/

	std::vector<double> log_x_vals;
	std::vector<double> log_h_vals;
	log_x_vals.reserve(N+2);
	log_h_vals.reserve(N+2);

	log_x_vals.push_back(-INFINITY);
	log_x_vals.push_back(log_L);
	log_x_vals.push_back(log_U);

	log_h_vals.push_back(log_prob_max);
	log_h_vals.push_back(log_p(log_L));
	log_h_vals.push_back(log_p(log_U));

	// We already have three (x, h(x)) pairs from above
	unsigned int iter = 3;

	while (!q.empty() && iter < N+2) {
		// Get the interval with the largest priority
		Interval int_top = q.top();
		q.pop();

		// Break the interval int_top into two pieces: left and right.
		double log_x_new = midpoint(int_top.get_log_x(), int_top.get_log_y());
		double log_h_new = log_p(log_x_new);

		// Add the midpoint to our list of knots
		log_x_vals.push_back(log_x_new);
		log_h_vals.push_back(log_h_new);
		iter++;

		// Add the interval which represents [int_top$log_x, log_x_new]
		Interval int_left(int_top.get_log_x(), log_x_new, int_top.get_log_h_x(),
			log_h_new, pw);
		q.push(int_left);

		// Add the interval which represents [log_x_new, int_top$log_x]
		Interval int_right(log_x_new, int_top.get_log_y(), log_h_new,
			int_top.get_log_h_y(), pw);
		q.push(int_right);

		Rcpp::checkUserInterrupt();
	}

	_log_x_vals.assign(log_x_vals.begin(), log_x_vals.end());
	_log_h_vals.assign(log_h_vals.begin(), log_h_vals.end());

	_knot_order = order(_log_x_vals, false);
	_log_x_vals.sort(false);
	_log_h_vals.sort(true);
}

inline const Rcpp::NumericVector& Stepdown::get_cum_probs() const
{
	return _cum_probs;
}

inline double Stepdown::get_norm_const() const
{
	return _norm_const;
}

inline const Rcpp::NumericVector& Stepdown::get_log_x_vals() const
{
	return _log_x_vals;
}

inline const Rcpp::NumericVector& Stepdown::get_log_h_vals() const
{
	return _log_h_vals;
}

inline const Rcpp::NumericVector& Stepdown::get_knot_order() const
{
	return _knot_order;
}

inline double Stepdown::quantile(double p, bool log) const
{
	// Recall that _cum_probs is a sorted vector
	// Use find_interval to locate the two cutpoints,
	// then do a linear interpolation between them.

	const Rcpp::NumericVector& zero = Rcpp::NumericVector::create(0);
	const Rcpp::NumericVector& cum_probs_ext = concat(zero, _cum_probs);

	unsigned int j1 = find_interval(p, cum_probs_ext);
	unsigned int j2 = j1 + 1;
	double cp1 = cum_probs_ext(j1);
	double cp2 = cum_probs_ext(j2);
	double log_x1 = _log_x_vals(j1);
	double log_x2 = _log_x_vals(j2);
	
	// Do the following computation, but be careful to keep x values on
	// the log-scale, since they may be extremely small.
	double log_ratio = std::log(p - cp1) - std::log(cp2 - cp1);
	double out = log_add2_exp(log_sub2_exp(log_x2, log_x1) + log_ratio, log_x1);

	return log ? out : std::exp(out);
}

inline double Stepdown::draw(bool log) const
{
	double u = R::runif(0, 1);
	return quantile(u, log);
}

inline Rcpp::NumericVector Stepdown::draw(unsigned int n, bool log) const
{
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = draw(log);
	}
	return out;
}

inline double Stepdown::density(double log_x, bool log, bool normalize) const
{
	// Get the idx such that log_h_vals[idx] <= log(x) < log_h_vals[idx+1]
	unsigned int idx = find_interval(log_x, _log_x_vals);
	double out = R_NegInf;

	if (idx < _N+2) {
		out = _log_h_vals(idx);
	}
	if (normalize) { out -= _norm_const; }
	return log ? out : exp(out);
}

inline double Stepdown::cdf(double log_x) const
{
	// Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]	
	unsigned int j1 = find_interval(log_x, _log_x_vals);
	unsigned int j2 = j1 + 1;
	double cp1 = _cum_probs(j1);
	double cp2 = _cum_probs(j2);
	double log_x1 = _log_x_vals(j1);
	double log_x2 = _log_x_vals(j2);
	
	// Do the following computation, but be careful to keep x values on
	// the log-scale, since they may be extremely small.
	double log_ratio = log_sub2_exp(log_x, log_x1) - log_sub2_exp(log_x2, log_x1) + std::log(cp2 - cp1);
	return exp(log_add2_exp(log(cp1), log_ratio));
}

// Note that pr_interval can be negative, even when x1 < x2, due to numerical
// accuracy. In this case, we will have to avoid taking the log.
inline double Stepdown::log_p(double log_u) const
{
	const std::pair<double,double>& endpoints = _w.roots(_w.log_c() + log_u);
	double log_pr = _g.pr_interval(endpoints.first, endpoints.second, true);

	// Be careful - it looks like we can get log(0) = NaN in C++
	return isnan(log_pr) ? R_NegInf : log_pr;
}

/*
* Update step function with new log_u and log_p values. This is currently
* inefficient for practical use. More sophisticated data structures would
* prevent the need to reallocate things.
*/
inline void Stepdown::add(double log_u)
{
	double log_h = log_p(log_u);

	_log_x_vals.push_back(log_u);
	_log_h_vals.push_back(log_h);
	_knot_order.push_back(_N + 3);

	// Need to subtract 1 from order results to get zero-based indices
	const Rcpp::IntegerVector& idx = order(_log_x_vals, false) - 1;

	_log_x_vals = _log_x_vals[idx];
	_log_h_vals = _log_h_vals[idx];
	_knot_order = _knot_order[idx];

	// Update cumulative probabilities
	_N++;
	update();
}

/*
* Compute normalizing constant, cumulative probabilities, and maximum jump
* based on current state of _log_x_vals and _log_h_vals. More sophisticated
* data structures would would make it possible to not reallocate things.
*/
inline void Stepdown::update()
{
	unsigned int N = _N;

	#ifdef SIMPLER_VERSION
		Rcpp::NumericVector areas(N+1);
		_cum_probs = Rcpp::NumericVector(N+1);
		_norm_const = 0;

		// This is a more readable version of the calculation that may be less precise
		for (unsigned int i = 0; i < N+1; i++) {
			double h_cur = exp(_log_h_vals(i));
			double h_next = exp(_log_h_vals(i+1));
			double x_cur = exp(_log_x_vals(i));
			double x_next = exp(_log_x_vals(i+1));
			areas(i) = h_cur * (x_next - x_cur);
			_norm_const += areas(i);
			_cum_probs(i) = _norm_const;
		}
	#else
		// This version tries to be more precise and do more work on the log-scale
		Rcpp::NumericVector log_areas(N+1);
		Rcpp::NumericVector log_cum_areas(N+1);

		log_areas(0) = _log_h_vals(0) + _log_x_vals(1);
		log_cum_areas(0) = log_areas(0);

		for (unsigned int l = 1; l < N+1; l++) {
			log_areas(l) = _log_h_vals(l) + log_sub2_exp(_log_x_vals(l+1), _log_x_vals(l));
			double arg1 = std::max(log_cum_areas(l-1), log_areas(l));
			double arg2 = std::min(log_cum_areas(l-1), log_areas(l));
			log_cum_areas(l) = log_add2_exp(arg1, arg2);
		}

		double log_normconst = log_cum_areas(N);
		const Rcpp::NumericVector& log_cum_probs = log_cum_areas - log_normconst;
		_norm_const = exp(log_normconst);
		_cum_probs = Rcpp::exp(log_cum_probs);
	#endif
}

}

#endif
