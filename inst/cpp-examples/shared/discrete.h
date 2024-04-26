#ifndef DISCRETE_H
#define DISCRETE_H

#include <Rcpp.h>

//' Quantile Function for Discrete Distributions with Finite Support
//' 
//' @param q A vector of quantiles to compute
//' @param cp Vector of cumulative probabilities of distribution
//' 
//' @return 0-based indices corresponding to the \code{q} quantiles
//' 
//' @details
//' Compute quantiles of a discrete distribution with cumulative probabilities
//' \code{cp(0)}, \code{cp(1)}, ..., \code{cp(k-1)}. Return indices in
//' \eqn{\{ 0, 1, \ldots, k-1 \}} which represent the
//' \code{q} quantiles. The caller can then identify the corresponding value
//' of the distribution (of which this function does not need to be aware).
//' 
//' Uses a bisection search so support relatively large \code{k}.
//' 
//' Note that both \code{q} and \code{cp} may be given on the original
//' probability scale or the log-scale, as long as they are both on the same
//' scale.
//' @examples
//' # Create a simple distribution and compute cumulative probabilities
//' x_seq = 1:10
//' k = length(x_seq)
//' p = rep(1/k, k)
//' cp = cumsum(p)
//' 
//' xx = sample(x = x_seq, size = 100000, replace = TRUE, prob = p)
//' 
//' Compare empirical qq quantiles with q_discrete
//' qq = seq(0, 1, length.out = 100)
//' idx_calc = q_discrete(qq, cp)
//' x_emp = quantile(xx, probs = qq)
//' plot(x_emp, x_seq[idx_calc + 1])
//' 
//' Repeat calculation on the log-scale
//' idx_calc_log = q_discrete(log(qq), log(cp))
//' plot(idx_calc, idx_calc_log)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector q_discrete(const Rcpp::NumericVector& q, const Rcpp::NumericVector& cp)
{
	unsigned int n = q.size();
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_discrete(q(i), cp);
	}

	return out;
}

unsigned int q_discrete(double q, const Rcpp::NumericVector& cp)
{
	unsigned int k = cp.size();
	if (q > cp(k-1)) {
		Rcpp::stop("q > max(cp)");
	}

	// Otherwise do a binary search
	// TBD: Could we use bisection functions here?
	unsigned int x_lo = 0;
	unsigned int x_hi = k-1;
	unsigned int x = (unsigned int)(floor((x_hi + x_lo) / 2.0));
	while (x_hi - x_lo > 1) {
		bool ind = (cp(x) >= q);
		x_lo = (1 - ind) * x + ind * x_lo;
		x_hi = ind * x + (1 - ind) * x_hi;
		x = (unsigned int)(floor((x_hi + x_lo) / 2.0));
	}

	if (cp(x_lo) >= q) {
		return x_lo;
	} else {
		return x_hi;
	}
}

#endif
