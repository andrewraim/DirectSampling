# ----- Setup for customized direct sampler -----
get_lognormal_weight = function(z, mu, sigma2)
{
	# log_c is the maximum value of the function log[w(x)]
	log_c = -(mu - sigma2 / 2)

	# Evaluate the function
	eval = function(x, log = FALSE) {
		n = length(x)
		out = rep(-Inf, n)
		idx = which(x < z)
		out[idx] = -log(z-x[idx]) - (log(z-x[idx]) - mu)^2 / (2*sigma2)
		if (log) { return(out) } else { return(exp(out)) }
	}

	# Return the roots of the equation w(x) = a, which is equivalent to
	# log(w(x)) = log(a). Roots are returned in increasing order.
	roots = function(log_a) {
		if (sigma2 < 2*(mu + log_a) && log_a <= log_c) {
			# Handle this case to prevent NaN issues due to numerical error
			x1 = z - exp(mu - sigma2)
			x2 = z - exp(mu - sigma2)
		} else {
			x1 = z - exp((mu - sigma2) + sqrt(sigma2 * (sigma2 - 2*(mu + log_a))))
			x2 = z - exp((mu - sigma2) - sqrt(sigma2 * (sigma2 - 2*(mu + log_a))))
		}

		# Edge case: If z is an integer and the larger root is numerically close to
		# z, make the root smaller by an epsilon. This helps to ensure that the
		# endpoint z is not part of the support of p(u).
		# if (z - floor(z) < 1e-10 && z - x2 < 1e-6) {
		#	x2 = z - 1e-6
		# }

		c(x1, x2)
	}

	ret = list(log_c = log_c, roots = roots, eval = eval)
	class(ret) = "weight"
	return(ret)
}
