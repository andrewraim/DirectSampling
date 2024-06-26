#include "direct_sampler_ar.h"
#include "Stepdown.h"

/*
* Note: This implementation is not complete. The calculation of M needs to
* be coded with the corrected expression.
*/
Rcpp::NumericVector direct_sampler_ar(unsigned int n, const WeightFunction& w,
	const BaseDistribution& g, double tol, unsigned int N_init,
	unsigned int max_rejections, const std::string& fill_method)
{
	// Use our Stepdown approximation to draw from p(u)
	Stepdown step(w, g, tol, N_init, fill_method);
	const Rcpp::NumericVector& v = step.draw(n);
	const Rcpp::NumericVector& log_u = Rcpp::log(v);

	Rcpp::NumericVector u(n);
	unsigned int rejections = 0;
	bool accept;

	// Because the step function is an upper bound for P(A_u), the constant M in
	// the acceptance ratio is always M = 1.
	double log_M = 0;
	
	for (unsigned int i = 0; i < n; i++) {
		accept = false;
		while (!accept && rejections < max_rejections) {
			double v = R::runif(0, 1);
			double u_proposal = step.draw_one();
			double log_p_val = step.log_p(log(u_proposal));
			double log_h_val = step.density(u_proposal, true, false);
			double log_ratio = log_p_val - log_h_val - log_M;
			if (log(v) < log_ratio) {
				// Accept u as a draw from p(u)
				u(i) = u_proposal;
				accept = true;
	 		} else {
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
		const std::pair<double,double>& endpoints = w.roots(w.log_c() + log(u(i)));
		out(i) = g.r_truncated(endpoints.first, endpoints.second);
		Rcpp::checkUserInterrupt();
	}

	return out;
}
