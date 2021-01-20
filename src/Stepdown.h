#ifndef STEPDOWN_H
#define STEPDOWN_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"

/*
* This class represents a step function approximation to the non-increasing
* Pr(A_u) for u in [0,1].
*/
class Stepdown {
public:
	Stepdown(const WeightFunction& w, const BaseDistribution& g,
		double tol, unsigned int N, const std::string& method);

	double get_norm_const() const;
	const Rcpp::NumericVector& get_cum_probs() const;
	const Rcpp::NumericVector& get_log_x_vals() const;
	const Rcpp::NumericVector& get_log_h_vals() const;

	void add(double log_u);
	double log_p(double log_u) const;
	double density(double x, bool take_log, bool normalize) const;
	double cdf(double x) const;
	double quantile(double p) const;
	double draw_one() const;
	Rcpp::NumericVector draw(unsigned int n) const;

private:
	void init_equal_steps(double log_L, double log_U, double log_prob_max,
		unsigned int N);
	void init_small_rects(double log_L, double log_U, double log_prob_max,
		unsigned int N);
	void update();

	const WeightFunction& _w;
	const BaseDistribution& _g;
	Rcpp::NumericVector _log_x_vals;
	Rcpp::NumericVector _log_h_vals;
	Rcpp::NumericVector _cum_probs;
	double _norm_const;
};

#endif
