#' Multivariate Normal Distribution
#' 
#' Density and drawing functions for the k-dimensional Multivariate Normal
#' distribution on, parameterized either by the covariance or precision
#' matrix.
#' 
#' @param n Number of draws to generate.
#' @param x An n-by-k matrix whose rows are points to evaluate.
#' @param mu k-dimensional vector representing the mean.
#' @param Sigma k-by-k covariance matrix.
#' @param Omega k-by-k precision matrix.
#' @param log If \code{TRUE}, return the log-density.
#' 
#' @name Multivariate Normal
NULL

#' @name Multivariate Normal
#' @export
rmvnorm = function(n, mu, Sigma)
{
	k = length(mu)
	stopifnot(all(dim(Sigma) == k))
	Z = matrix(rnorm(n*k), k, n)
	A = t(chol(Sigma))
	A %*% Z + mu
}

#' @name Multivariate Normal
#' @export
rmvnorm_prec = function(n, mu, Omega)
{
	k = length(mu)
	stopifnot(all(dim(Omega) == k))
	Z = matrix(rnorm(n*k), k, n)

	# Note that Ainv %*% t(Ainv) is the Cholesky decomposition of the covariance
	# matrix solve(Omega)
	A = chol(Omega)
	Ainv = backsolve(A, diag(1,k,k))
	Ainv %*% Z + mu
}

#' @name Multivariate Normal
#' @export
dmvnorm = function(x, mu, Sigma, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(all(dim(Sigma) == k))
	stopifnot(ncol(x) == k)

	# Use QR decomposition to efficiently compute determinant and inverse
	xc = x - t(mu) %x% matrix(1,n,1)
	qr_out = qr(Sigma)
	Omega_xc = qr.solve(qr_out, t(xc))
	logdetA = sum(log(abs(diag(qr.R(qr_out)))))
	logf = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(logf) } else { return(exp(logf))}
}

#' @name Multivariate Normal
#' @export
dmvnorm_prec = function(x, mu, Omega, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(all(dim(Omega) == k))
	stopifnot(ncol(x) == k)

	xc = x - t(mu) %x% matrix(1,n,1)
	Omega_xc = Omega %*% t(xc)
	logdetA = -as.numeric(determinant(Omega)$modulus)
	logf = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(logf) } else { return(exp(logf))}
}
