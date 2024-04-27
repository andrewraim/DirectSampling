#ifndef DIRECT_SAMPLER_AR_H
#define DIRECT_SAMPLER_AR_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"
#include "Stepdown.h"

namespace DirectSampling {

Rcpp::List direct_sampler_ar(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N,
	unsigned int max_rejections, const std::string& fill_method,
	double priority_weight = 0.5)
{
	// Use our Stepdown approximation to draw from p(u)
	Stepdown step(w, g, tol, N, fill_method, priority_weight);

	Rcpp::NumericVector log_u(n);
	unsigned int rejections = 0;
	bool accept;

	// Because the step function is an upper bound for P(A_u), the constant M in
	// the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n; i++) {
		accept = false;
		while (!accept && rejections < max_rejections) {
			double v = R::runif(0, 1);
			double log_u_proposal = step.draw_one(true);
			double log_p_val = step.log_p(log_u_proposal);
			double log_h_val = step.density(log_u_proposal, true, false);
			double log_ratio = log_p_val - log_h_val - log_M;

			if (log(v) < log_ratio) {
				// Accept u as a draw from p(u)
				log_u(i) = log_u_proposal;
				accept = true;
	 		} else {
	 			// Reject u and add it to knots
	 			step.add(log_u_proposal);
				rejections++;
			}
	 		Rcpp::checkUserInterrupt();
		}
	}

	if (rejections == max_rejections) {
		char msg[64];
		sprintf(msg, "Reached maximum number of rejections: %d\n", max_rejections);
		Rcpp::stop(msg);
	}

	// Draw from g(x | u) for each value of log(u)
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log_u(i));
		out(i) = g.r_truncated(endpoints.first, endpoints.second);
		Rcpp::checkUserInterrupt();
	}

	return Rcpp::List::create(
		Rcpp::Named("x") = out,
		Rcpp::Named("rejections") = rejections
	);
}

}

#endif
