# Basic direct sampler, without any of the customizations we proposed
#' @export
direct_sampler_basic = function(n, w_obj, g_obj, N = 100, max_iter = 10000)
{
	# Maximum weight c obtained from the mode of the LogNormal distn.
	log_c = w_obj$log_c

	# Set up approximation for the marginal distn of u
	knots = 0:N / N
	qq = numeric(N+1)
	for (k in 1:length(qq)) {
		# Pr(A_{k/N})
		u = knots[k]
		endpoints = w_obj$roots(log_c + log(u))
		qq[k] =  g_obj$pr_interval(endpoints[1], endpoints[2])
	}

	# Draw from approximation to marginal distn of u
	k = sample(0:N, size = n, replace = TRUE, prob = qq)
	u = rbeta(n, k+1, N-k+1)

	x = rep(NA, n)
	tries = numeric(n)

	for (i in 1:n) {
		while (is.na(x[i]) && tries[i] < max_iter) {
			tries[i] = tries[i] + 1
			x_proposal = g_obj$r_truncated(1)
			accept = w_obj$eval(x_proposal, log = TRUE) > log_c + log(u[i])
			if (accept) { x[i] = x_proposal }
		}
	}

	list(x = x, tries = tries)
}

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
