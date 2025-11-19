#ifndef DIRECT_SAMPLING_REJECTION_H
#define DIRECT_SAMPLING_REJECTION_H

#include <Rcpp.h>
#include "WeightFunction.h"
#include "BaseDistribution.h"
#include "Stepdown.h"
#include "rejection-args.h"

namespace DirectSampling {

inline Rcpp::List rejection(Stepdown step, unsigned int n,
	const rejection_args& args = rejection_args())
{
	unsigned int max_rejects = args.max_rejects;
	bool adapt = args.adapt;
	unsigned int report = args.report;
	fntl::error_action action = args.action;

	Rcpp::NumericVector log_u(n);
	const WeightFunction& w = step.get_weight();
	const BaseDistribution& g = step.get_base();

	std::vector<unsigned int> rejects;
	unsigned int N_rejects = 0;

	// Because the step function is an upper bound for P(A_u), the constant M in
	// the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n; i++) {
		bool accept = false;
		rejects.push_back(0L);
		
		while (!accept && N_rejects < max_rejects) {
			double v = R::runif(0, 1);
			double log_u_proposal = step.draw(true);
			double log_p_val = step.log_p(log_u_proposal);
			double log_h_val = step.density(log_u_proposal, true, false);
			double log_ratio = log_p_val - log_h_val - log_M;

			if (log(v) < log_ratio) {
				// Accept u as a draw from p(u)
				log_u(i) = log_u_proposal;
				accept = true;
			} else if (adapt) {
	 			// Reject u and add it to knots
	 			step.add(log_u_proposal);
				N_rejects++;
				rejects[i]++;
			} else {
				N_rejects++;
				rejects[i]++;
			}
			
			// Report progress after `report` candidates
			unsigned int N_accepts = i + accept;
			if ((N_rejects + N_accepts) % report == 0) {
				Rprintf("%s - After %d candidates, %d accepts and %d rejects\n",
					timestamp().c_str(), N_accepts + N_rejects, N_accepts, N_rejects);
			}

			if ((N_rejects + N_accepts) % 10000 == 0) {
		 		Rcpp::checkUserInterrupt();
			}
		}
	}

	if (N_rejects > max_rejects) {
		switch(action) {
			case fntl::error_action::STOP:
				Rcpp::stop("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case fntl::error_action::WARNING:
				Rcpp::warning("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case fntl::error_action::MESSAGE:
				Rprintf("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			default:
				break;
		}
	}

	// Draw from g(x | u) for each value of log(u)
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log_u(i));
		out(i) = g.r_truncated(endpoints.first, endpoints.second);

		if (i % 10000 == 0) {
		 	Rcpp::checkUserInterrupt();
		}
	}

	return Rcpp::List::create(
		Rcpp::Named("x") = out,
		Rcpp::Named("rejects") = rejects
	);
}

}

#endif
