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