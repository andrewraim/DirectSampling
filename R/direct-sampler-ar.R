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
#' @name Direct Sampler AR
NULL

#' @name Direct Sampler AR
#' @export
direct_sampler_ar = function(n, w, g, tol = 1e-8, N = 10,
	max_rejections, fill_method = "small_rects")
{
	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Use our Stepdown approximation to draw from p(u)
	u = numeric(n)
	step = Stepdown$new(w, g, tol, N = N, fill_method)
	rejections = 0

	# Because the step function is an upper bound for P(A_u), the constant M in
	# the acceptance ratio is always M = 1.
	log_M = 0

	for (i in 1:n) {
		accept = FALSE
		while (!accept && rejections < max_rejections) {
			v = runif(1)
			u_proposal = step$r(1)
			log_p_val = step$get_log_p(log(u_proposal))
			log_h_val = step$d(u_proposal, log = TRUE, normalize = FALSE)
			log_ratio = log_p_val - log_h_val - log_M
			if (log(v) < log_ratio) {
				# Accept u as a draw from p(u)
				u[i] = u_proposal
				accept = TRUE
			} else {
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

direct_sampler_ar_adapt = function(n, w, g, tol = 1e-8, N = 10,
	max_rejections, fill_method = "small_rects")
{
	stop("This function is under construction!")

	stopifnot(class(w) == "weight")
	stopifnot(class(g) == "base")

	# Use our Stepdown approximation to draw from p(u)
	u = numeric(n)
	step = Stepdown$new(w, g, tol, N = N, fill_method)
	rejections = 0

	# Because the step function is an upper bound for P(A_u), the constant M in
	# the acceptance ratio is always M = 1.
	log_M = 0

	for (i in 1:n) {
		accept = FALSE
		while (!accept && rejections < max_rejections) {
			v = runif(1)
			u_proposal = step$r(1)
			log_p_val = step$get_log_p(log(u_proposal))
			log_h_val = step$d(u_proposal, log = TRUE, normalize = FALSE)
			log_ratio = log_p_val - log_h_val - log_M
			if (log(v) < log_ratio) {
				# Accept u as a draw from p(u)
				u[i] = u_proposal
				accept = TRUE
	 		} else {
				# TBD: Adapt the step function.
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