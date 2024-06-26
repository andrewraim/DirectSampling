#' Direct sampler with accept-reject
#'
#' @param n Number of draws to generate.
#' @param w A weight function object.
#' @param g A base distribution object.
#' @param tol Tolerance for step function approximation in customized sampler.
#' @param N Number of knots to use in approximation for \eqn{p(u)}. Specifically,
#' N+1 points will be used.
#' @param max_rejections Maximum number of rejections to allow. If this number
#' is reached, an exception will be thrown.
#' @param fill_method Knot selection method for customized sampler. See
#' \code{\link{Stepdown}}.
#' @param priority_weight TBD
#' @param verbose Return a list with the sample \code{x} and the number of
#' rejections \code{rejections}.
#'
#' @return A vector of length \code{n} which represents the sample.
#' 
#' @details
#' Direct sampling with accept-reject algorithm can be used to take exact draws
#' from the target distribution, with the added cost that same candidate draws
#' will be rejected.
#' 
#' The function \code{direct_sampler_ar} does not adapt the underlying step
#' function using rejections; here, the number of knots \code{N} is kept
#' constant throughout sampling.
#' 
#' @name direct-sampler-ar
#' @export
direct_sampler_ar = function(n, w, g, tol = 1e-8, N = 10,
	max_rejections, fill_method = "small_rects", priority_weight = 1/2,
	verbose = FALSE)
{
	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Use our Stepdown approximation to draw from p(u)
	log_u = numeric(n)
	step = Stepdown$new(w, g, tol, N, fill_method, priority_weight)
	rejections = 0

	# Because the step function is an upper bound for P(A_u), the constant M in
	# the acceptance ratio is always M = 1.
	log_M = 0

	for (i in 1:n) {
		accept = FALSE
		while (!accept && rejections < max_rejections) {
			v = runif(1)
			log_u_proposal = step$r(1, log = TRUE)
			log_p_val = step$get_log_p(log_u_proposal)
			log_h_val = step$d(log_u_proposal, log = TRUE, normalize = FALSE)
			log_ratio = log_p_val - log_h_val - log_M

			if (log(v) < log_ratio) {
				# Accept u as a draw from p(u)
				log_u[i] = log_u_proposal
				accept = TRUE
			} else {
				step$add(log_u_proposal)
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
		endpoints = w$roots(w$log_c + log_u[i])
		x[i] = g$r_truncated(1, endpoints[1], endpoints[2])
	}

	if (verbose) {
		out = list(x = x, rejections = rejections)
	} else {
		out = x
	}
	return(out)
}

