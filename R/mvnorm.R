#' Multivariate Normal Distribution
#' 
#' Density and drawing functions for Multivariate Normal distribution on
#' \eqn{\mathbb{R}^k}, parameterized either by the covariance or precision
#' matrix.
#' 
#' @param n Number of draws to generate.
#' @param x An n-by-k matrix whose rows are points to evaluate.
#' @param mu k-dimensional vector representing the mean.
#' @param Sigma k-by-k covariance matrix.
#' @param Omega k-by-k precision matrix.
#' @param log If \code{TRUE}, return the log-density.
#' 
#' @name MVN
NULL

#' @name MVN
#' @export
rmvnorm = function(n, mu, Sigma)
{
	k = length(mu)
	stopifnot(k == nrow(Sigma) && k == ncol(Sigma))
	Z = matrix(rnorm(n*k), k, n)
	A = t(chol(Sigma))
	A %*% Z + mu
}

#' @name MVN
#' @export
rmvnorm_prec = function(n, mu, Omega)
{
	k = length(mu)
	stopifnot(k == nrow(Omega) && k == ncol(Omega))
	Z = matrix(rnorm(n*k), k, n)

	# Note that Ainv %*% t(Ainv) is the Cholesky decomposition of the covariance
	# matrix solve(Omega)
	A = chol(Omega)
	Ainv = backsolve(A, diag(1,k,k))
	Ainv %*% Z + mu
}

#' @name MVN
#' @export
dmvnorm = function(x, mu, Sigma, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(k == nrow(Sigma) && k == ncol(Sigma) && k == ncol(x))

	# Use QR decomposition to efficiently compute determinant and inverse
	xc = x - t(mu) %x% matrix(1,n,1)
	qr_out = qr(Sigma)
	Omega_xc = qr.solve(qr_out, t(xc))
	logdetA = sum(log(abs(diag(qr.R(qr_out)))))
	logf = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(logf) } else { return(exp(logf))}
}

#' @name MVN
#' @export
dmvnorm_prec = function(x, mu, Omega, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(k == nrow(Sigma) && k == ncol(Sigma) && k == ncol(x))

	xc = x - t(mu) %x% matrix(1,n,1)
	Omega_xc = Omega %*% t(xc)
	logdetA = -as.numeric(determinant(Omega)$modulus)
	logf = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(logf) } else { return(exp(logf))}
}

