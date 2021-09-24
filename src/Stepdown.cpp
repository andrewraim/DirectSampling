#include "Stepdown.h"
#include <math.h>
#include <queue>
#include "util.h"
#include "functionals.h"
#include "bisection.h"
#include "find_interval.h"

/*
* A predicate which is used to locate log_L
*/
class Predicate_logL : public Predicate
{
public:
	Predicate_logL(const Stepdown& step, double log_prob_max)
		: _step(step), _log_prob_max(log_prob_max)
	{
	}
	bool operator()(double log_u) const {
		double log_prob = _step.log_p(log_u);
		return log_prob < _log_prob_max;
	}
private:
	const Stepdown& _step;
	double _log_prob_max;
};

/*
* A predicate which is used to locate log_U
*/
class Predicate_logU : public Predicate
{
public:
	Predicate_logU(const Stepdown& step, double log_prob_max)
		: _step(step), _log_prob_max(log_prob_max)
	{
	}
	bool operator()(double log_u) const {
		double log_prob = _step.log_p(log_u);
		// return std::isinf(log_prob);
		return log_prob < log(1e-5) + _log_prob_max;
	}
private:
	const Stepdown& _step;
	double _log_prob_max;
};

/*
* This class is used in the heap implementation of Stepdown::init_small_rects.
* It represents an interval [x,y] with function values h(x) and h(y).
* We save width = y-x and height = h(x) - h(y).
*/
class Interval
{
public:
	Interval(double log_x, double log_y, double log_h_x, double log_h_y, double pw = 0.5)
		: _log_x(log_x), _log_y(log_y), _log_h_x(log_h_x), _log_h_y(log_h_y), _pw(pw)
	{
		if (_pw < 0 || _pw > 1) {
			Rcpp::stop("pw must be in [0,1]");
		}
	}
	Interval(const Interval& rhs) {
		*this = rhs;
	}
	double width() const {
		return exp(_log_y) - exp(_log_x);
	}
	double height() const {
		return exp(_log_h_x) - exp(_log_h_y);
	}
	double priority() const {
		double width = exp(_log_y) - exp(_log_x);
		double height = exp(_log_h_x) - exp(_log_h_y);
		return _pw*log(height) + (1-_pw)*log(width);
	}
	const Interval& operator=(const Interval& rhs) {
		_log_x = rhs._log_x;
		_log_y = rhs._log_y;
		_log_h_x = rhs._log_h_x;
		_log_h_y = rhs._log_h_y;
		_pw = rhs._pw;
		return *this;
	}
	bool operator<(const Interval& rhs) const {
		return priority() < rhs.priority();
	}
	bool operator>(const Interval& rhs) const {
		return priority() > rhs.priority();
	}
	double get_log_x() const {
		return _log_x;
	}
	double get_log_y() const {
		return _log_y;
	}
	double get_log_h_x() const {
		return _log_h_x;
	}
	double get_log_h_y() const {
		return _log_h_y;
	}
	void print(bool log_scale = false) const {
		if (log_scale) {
			Rprintf("log_x: %g\n", _log_x);
			Rprintf("log_y: %g\n", _log_y);
			Rprintf("log_h_x: %g\n", _log_h_x);
			Rprintf("log_h_y: %g\n", _log_h_y);
			Rprintf("width: %g\n", width());
			Rprintf("height: %g\n", height());
			Rprintf("priority: %g\n", priority());
		} else {
			Rprintf("x: %g\n", exp(_log_x));
			Rprintf("y: %g\n", exp(_log_y));
			Rprintf("h_x: %g\n", exp(_log_h_x));
			Rprintf("h_y: %g\n", exp(_log_h_y));
			Rprintf("width: %g\n", width());
			Rprintf("height: %g\n", height());
			Rprintf("priority: %g\n", priority());
		}
	}
private:
	double _log_x;
	double _log_y;
	double _log_h_x;
	double _log_h_y;
	double _pw;
};

void print(std::priority_queue<Interval> q)
{
	while(!q.empty()) {
		Interval el = q.top();
		q.pop();
		el.print();
	}
}

Stepdown::Stepdown(const WeightFunction& w, const BaseDistribution& g,
	double tol, unsigned int N, const std::string& method, double priority_weight)
	: _w(w), _g(g), _tol(tol), _N(N), _log_x_vals(), _log_h_vals(),
	  _cum_probs(), _priority_weight(priority_weight)
{
	double log_prob;
	double delta;
	std::pair<double,double> endpoints;
	Midpoint midpoint;
	IntervalLength dist;
	
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
	Predicate_logL predL(*this, log_prob_max);
	delta = tol * (log_L_hi - log_L_lo);
	double log_L = bisection(log_L_lo, log_L_hi, predL, midpoint, dist, delta);

	// Do a bisection search to find U, the smallest point where P(A_U) = 0.
	Predicate_logU predU(*this, log_prob_max);
	delta = tol * (0 - log_L);
	double log_U = bisection(log_L, 0, predU, midpoint, dist, delta);

	// Now fill in points between L and U
	if (method == "equal_steps") {
		init_equal_steps(log_L, log_U, log_prob_max);
	} else if (method == "small_rects") {
		init_small_rects(log_L, log_U, log_prob_max);
	} else {
		char msg[256];
		sprintf(msg,
          "Unknown method: %s. Currently support equal_steps and small_rects\n",
          method.c_str());
		Rcpp::stop(msg);
	}

	update();
}

void Stepdown::init_equal_steps(double log_L, double log_U, double log_prob_max)
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
}

void Stepdown::init_small_rects(double log_L, double log_U, double log_prob_max)
{
	unsigned int N = _N;
	double pw = _priority_weight;
	Midpoint midpoint;

	// This queue should be in max-heap order by area
	std::priority_queue<Interval> q;
	q.push(Interval(log_L, log_U, log_p(log_L), log_p(log_U), pw));

	// Try to be efficient by preallocating log_x_vals and
	// log_h_vals to the (maximum) size needed, then copying these temporary
	// structures to _log_x_vals and _log_h_vals afterward. (It looks like
	// Rcpp NumericVectors cannot have their allocation controlled in this way).
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
		
		// Rprintf("-----------------\n");
		// Rprintf("Popped int_top:\n");
		// int_top.print();
		// Rprintf("\n");

		// Break the interval int_top into two pieces: left and right.
		double log_x_new = log(midpoint(exp(int_top.get_log_x()), exp(int_top.get_log_y())));
		double log_h_new = log_p(log_x_new);

		// Add the midpoint to our list of knots
		log_x_vals.push_back(log_x_new);
		log_h_vals.push_back(log_h_new);
		iter++;
		// if (log_h_new > log(1e-5) + log_prob_max) {
			// Experimental: Only add x_new if h_new > _tol*prob_max
			// Trying to keep inconsequential knot points out of final set
			// log_x_vals.push_back(log_x_new);
			// log_h_vals.push_back(log_h_new);
			// iter++;
		// }

		// Add the interval which represents [int_top$log_x, log_x_new]
		Interval int_left(int_top.get_log_x(), log_x_new, int_top.get_log_h_x(),
			log_h_new, pw);
		q.push(int_left);

		// Add the interval which represents [log_x_new, int_top$log_x]
		Interval int_right(log_x_new, int_top.get_log_y(), log_h_new,
			int_top.get_log_h_y(), pw);
		q.push(int_right);

		// Rprintf("Pushed int_left:\n");
		// int_left.print();
		// Rprintf("\n");
		// Rprintf("Pushed int_right:\n");
		// int_right.print();

		Rcpp::checkUserInterrupt();
	}

	_log_x_vals.assign(log_x_vals.begin(), log_x_vals.end());
	_log_h_vals.assign(log_h_vals.begin(), log_h_vals.end());
	_log_x_vals.sort(false);
	_log_h_vals.sort(true);

	// Rprintf("After init:\n");
	// Rcpp::NumericVector x_vals = Rcpp::exp(_log_x_vals);
	// Rcpp::NumericVector h_vals = Rcpp::exp(_log_h_vals);
	// Rcpp::print(x_vals);
	// Rcpp::print(h_vals);
}

const Rcpp::NumericVector& Stepdown::get_cum_probs() const
{
	return _cum_probs;
}

double Stepdown::get_norm_const() const
{
	return _norm_const;
}

const Rcpp::NumericVector& Stepdown::get_log_x_vals() const
{
	return _log_x_vals;
}

const Rcpp::NumericVector& Stepdown::get_log_h_vals() const
{
	return _log_h_vals;
}

double Stepdown::quantile(double p) const
{
	// Recall that _cum_probs is a sorted vector
	// Use find_interval to locate the two cutpoints,
	// then do a linear interpolation between them.

	const Rcpp::NumericVector& zero = Rcpp::NumericVector::create(0);
	const Rcpp::NumericVector& cum_probs_ext = concat(zero, _cum_probs);

	unsigned int j1 = find_interval(p, cum_probs_ext);
	unsigned int j2 = j1 + 1;
	double x1 = exp(_log_x_vals(j1));
	double x2 = exp(_log_x_vals(j2));
	double cp1 = cum_probs_ext(j1);
	double cp2 = cum_probs_ext(j2);
	return x1 + (x2 - x1) * (p - cp1) / (cp2 - cp1);
}

double Stepdown::draw_one() const
{
	double u = R::runif(0, 1);
	return quantile(u);
}

Rcpp::NumericVector Stepdown::draw(unsigned int n) const
{
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = draw_one();
	}
	return out;
}

double Stepdown::density(double x, bool take_log, bool normalize) const
{
	// Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]
	unsigned int idx = find_interval(log(x), _log_x_vals);
	double out = _log_h_vals(idx);
	if (normalize) { out -= _norm_const; }
	if (take_log) { return out; } else { return exp(out); }
}

double Stepdown::cdf(double x) const
{
	// Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]	
	unsigned int j1 = find_interval(log(x), _log_x_vals);
	unsigned int j2 = j1 + 1;
	double x1 = exp(_log_x_vals(j1));
	double x2 = exp(_log_x_vals(j2));
	double cp1 = _cum_probs(j1);
	double cp2 = _cum_probs(j2);
	return cp1 + (x - x1) / (x2 - x1) * (cp2 - cp1);
}

double Stepdown::log_p(double log_u) const
{
	const std::pair<double,double>& endpoints = _w.roots(_w.log_c() + log_u);
	return log(_g.pr_interval(endpoints.first, endpoints.second));
}

/*
* Update step function with new log_u and log_p values. This is currently
* inefficient for practical use. More sophisticated data structures would
* prevent the need to reallocate things.
*/
void Stepdown::add(double log_u)
{
	double log_h = log_p(log_u);

	// This vector is sorted in nondecreasing order	
	_log_x_vals.push_back(log_u);
	_log_x_vals.sort(false);

	// This vector is sorted in nonincreasing order	
	_log_h_vals.push_back(log_h);
	_log_h_vals.sort(true);

	update();
}

/*
* Compute normalizing constant, cumulative probabilities, and maximum jump
* based on current state of _log_x_vals and _log_h_vals. More sophisticated
* data structures would would make it possible to not reallocate things.
*/
void Stepdown::update()
{
	unsigned int N = _log_x_vals.size() - 2;
	Rcpp::NumericVector areas(N+1);

	_cum_probs = Rcpp::NumericVector(N+1);
	_norm_const = 0;

	#ifdef SIMPLER_VERSION
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
		for (unsigned int i = 0; i < N+1; i++) {
			double log_h_cur = _log_h_vals(i);
			// double log_h_next = _log_h_vals(i+1);
			double log_x_cur = _log_x_vals(i);
			double log_x_next = _log_x_vals(i+1);
			areas(i) = exp(log_x_next + log_h_cur) - exp(log_x_cur + log_h_cur);
			_norm_const += areas(i);
			_cum_probs(i) = _norm_const;
		}
	#endif

	_cum_probs = _cum_probs / _norm_const;

	if (Rcpp::is_true(Rcpp::any(Rcpp::is_na(_cum_probs)))) {
		Rcpp::warning("NA values found in cum_probs. Approximation may have failed");
	}
}
