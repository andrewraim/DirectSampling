#' Inverse Gamma Distribution
#' 
#' Distribution functions for Inverse Gamma.Simplified versions of code from
#' the invgamma package.

#' 
#' @param n Number of draws to generate
#' @param x A vector of n points to evaluate.
#' @param a shape parameter.
#' @param b rate parameter.
#' @param log If \code{TRUE} return the log-density.
#' 
#' @details
#' Here we assume Inverse Gamma distribution with density
#' \deqn{
#' f(x) = \frac{b^a}{\Gamma(a)} x^{-(a+1)} e^{-b/x}.
#' }
#' 
#' @name Inverse Gamma
NULL

#' @name Inverse Gamma
#' @export
rinvgamma = function(n, a, b)
{
	1 / rgamma(n, a, b)
}

#' @name Inverse Gamma
#' @export
dinvgamma = function(x, a, b, log = FALSE)
{
	logf = dgamma(1/x, a, b, log = TRUE) - 2 * log(x)
	if (log) { return(logf) } else { return(exp(logf))}
}

#' @name Inverse Gamma
#' @export
pinvgamma = function(x, a, b) {
	pgamma(1/x, a, b, lower.tail = FALSE)
}

#' @name Inverse Gamma
#' @export
qinvgamma = function(p, a, b) {
	1 / qgamma(1-p, a, b)
}
