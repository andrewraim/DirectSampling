#' Dirichlet Distribution
#' 
#' Density and drawing function for Dirichlet distribution on the \eqn{k}
#' dimensional probability simplex.
#' 
#' @param n Number of draws to generate
#' @param x An n-by-k matrix whose rows are points to evaluate.
#' @param alpha k-dimensional vector of intensity parameters.
#' @param log If \code{TRUE} return the log-density.
#' 
#' @details
#' Here we assume Dirichlet distribution with density
#' \deqn{
#' f(x) = \frac{x_1^{\alpha_1-1} \cdots x_k^{\alpha_k-1}}{B(\bm{\alpha})}.
#' }
#' 
#' @name Dirichlet
NULL

#' @name Dirichlet
#' @export
rdirichlet = function(n, alpha)
{
	k = length(alpha)
	x = matrix(NA, n, k)
	for (j in 1:k) {
		x[,j] = rchisq(n, 2*alpha[j])
	}
	S = rowSums(x) %*% t(rep(1,k))
	return(x / S)
}

#' @name Dirichlet
#' @export
ddirichlet = function(x, alpha, log = FALSE)
{
	n = nrow(x)
	k = ncol(x)
	stopifnot(k == length(alpha))
	logf = rep(lgamma(sum(alpha)) - sum(lgamma(alpha)), n)

	for (j in 1:k) {
		logf = logf + (alpha[j]-1)*x[,j]
	}

	if (log) { return(logf) } else { return(exp(logf))}
}
