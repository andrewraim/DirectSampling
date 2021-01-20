rinvgamma = function(n, a, b)
{
	1 / rgamma(n, a, b)
}

dinvgamma = function(x, a, b, log = FALSE)
{
	logf = dgamma(1/x, a, b, log = TRUE) - 2 * log(x)
	if (log) { return(logf) } else { return(exp(logf))}
}