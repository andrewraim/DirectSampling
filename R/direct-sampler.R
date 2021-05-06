#' Direct Sampler
#' 
#' Basic and customized versions of the direct sampler.
#' 
#' @param n Number of draws to generate.
#' @param w A weight function object.
#' @param g A base distribution object.
#' @param N Number of knots to use in approximation for \eqn{p(u)}. Specifically,
#' N+1 points will be used.
#' @param max_iter Maximum number of accept/reject attempts to try for each draw.
#' @param tol Tolerance for step function approximation in customized sampler.
#' @param fill_method Knot selection method for customized sampler. See
#' \code{\link{Stepdown}}.
#' 
#' @name Direct Sampler
NULL

#' @name Direct Sampler
#' @export
direct_sampler_basic = function(n, w, g, N = 100, max_iter = 10000)
{
	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Maximum weight c obtained from the mode of the LogNormal distn.
	log_c = w$log_c

	# Set up approximation for the marginal distn of u
	knots = 0:N / N
	qq = numeric(N+1)
	for (k in 1:length(qq)) {
		# Pr(A_{k/N})
		u = knots[k]
		endpoints = w$roots(log_c + log(u))
		qq[k] =  g$pr_interval(endpoints[1], endpoints[2])
	}

	# Draw from approximation to marginal distn of u
	k = sample(0:N, size = n, replace = TRUE, prob = qq)
	u = rbeta(n, k+1, N-k+1)

	x = rep(NA, n)
	tries = numeric(n)

	for (i in 1:n) {
		while (is.na(x[i]) && tries[i] < max_iter) {
			tries[i] = tries[i] + 1
			x_proposal = g$r_truncated(1)
			accept = w$eval(x_proposal, log = TRUE) > log_c + log(u[i])
			if (accept) { x[i] = x_proposal }
		}
	}

	list(x = x, tries = tries)
}

#' @name Direct Sampler
#' @export
direct_sampler = function(n, w, g, tol = 1e-8, N = 100,
	fill_method = "small_rects")
{
	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Use our Stepdown approximation to draw from p(u)
	step = Stepdown$new(w, g, tol, N, fill_method)
	u = step$r(n)

	# Draw from g(x | u) for each value of u
	x = rep(NA, n)
	for (i in 1:n) {
		endpoints = w$roots(w$log_c + log(u[i]))
		x[i] = g$r_truncated(1, endpoints[1], endpoints[2])
	}

	return(x)
}
