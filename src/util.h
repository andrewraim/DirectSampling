#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

/*
*  The following macros help us to use the properties:
*  log(x + y) = log(x) + log(1 + y/x)
*             = log(x) + log1p(exp(log(y) - log(x))),
*  log(x - y) = log(x) + log(1 - y/x)
*             = log(x) + log1p(-exp(log(y) - log(x))),
*  with log(x) and log(y) given as inputs, i.e. on the log-scale. When x and y
*  are of very different magnitudes, this is more stable when x is taken to be
*  the larger of the inputs. The most extreme case is when one of the inputs
*  might be -inf; in this case that input should be the second one.
*
*  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
*/
#define logadd(logx, logy) (logx + log1p(exp(logy - logx)))
#define logsub(logx, logy) (logx + log1p(-exp(logy - logx)))

/*
* These didn't exist in Rcpp Sugar (as far as I could tell), so we'll make our
* own quick versions.
*/
Rcpp::NumericVector concat(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);
Rcpp::IntegerVector order(const Rcpp::NumericVector& x, bool decrease);
Rcpp::NumericVector cumsum(const Rcpp::NumericVector& x);

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
Rcpp::IntegerVector q_discrete(const Rcpp::NumericVector& q, const Rcpp::NumericVector& cp);

unsigned int q_discrete(double q, const Rcpp::NumericVector& cp);

#endif
