direct_sampler_ar = function(n, w, g, tol = 1e-8, N_init = 10,
	max_rejections, fill_method = "small_rects")
{
	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Use our Stepdown approximation to draw from p(u)
	u = numeric(n)
	step = Stepdown$new(w, g, tol, N = N_init, fill_method)
	rejections = 0

	for (i in 1:n) {
		accept = FALSE
		while (!accept && rejections < max_rejections) {
			v = runif(1)
			u_proposal = step$r(1)
			log_p_val = step$get_log_p(log(u_proposal))
			log_h_val = step$d(u_proposal, log = TRUE, normalize = FALSE)
			M = step$get_max_jump()
			log_ratio = log_p_val - log_h_val - log(M)
			if (log(v) < log_ratio) {
				# Accept u as a draw from p(u)
				u[i] = u_proposal
				accept = TRUE
	 		} else {
				# Adapt the step function
				# printf("Adding %f to step function\n", log(u))
				step$add(log(u_proposal))
				rejections = rejections + 1
			}
		}
	}

	if (rejections == max_rejections) {
		msg = sprintf("Reached maximum number of rejections: %d\n", max_rejections)
		stop(msg)
	}

	# Draw from g(x | u) for each value of log(u)
	x = rep(NA, n)
	for (i in 1:n) {
		endpoints = w$roots(w$log_c + log(u[i]))
		x[i] = g$r_truncated(1, endpoints[1], endpoints[2])
	}

	return(x)
}