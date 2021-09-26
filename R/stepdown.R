#' @title Stepdown class
#' 
#' @param w An object representing a weight function (see details).
#' @param g An object representing a base distribution (see details).
#' @param tol A small positive number used in search for \eqn{u_L} and \eqn{u_H}.
#' @param N Number of knot points will be N+1.
#' @param method Either \code{equal_steps} or \code{small_rects} (see details).
#' @param log_u A value of u given at the log-scale.
#' @param p A probability to evaluate.
#' @param n Number of draws to generate.
#' @param x A vector of points to evaluate.
#' @param log If \code{TRUE}, return value on the log-scale.
#' @param normalize If \code{TRUE}, normalize the result to be a density value.
#' 
#' @details
#' This R6 class represents a step function approximation to the non-increasing
#' function \eqn{\rm{Pr}(A_u)} for \eqn{u \in [0,1]}.
#' 
#' The object \code{w} represents a unimodal weight function. It is
#' expected to contain several members:
#' \itemize{
#' \item \code{log_c}: the logarithm of the value c, which is the mode of the
#' weight function.
#' \item \code{roots(log_a)}: return the roots of the equation
#' \eqn{\log w(x) - \rm{log_a} = 0}.
#' \item \code{eval(x, log = TRUE)}: evaluate the function. Return the value
#' on the log-scale if \code{log = TRUE}.
#' }
#' 
#' The object \code{g} represents a base distribution. It is
#' expected to contain two member functions:
#' \itemize{
#' \item \code{pr_interval(x1, x2)}: Return \eqn{\rm{Pr}(x1 < X < x2)} under
#' distribution g.
#' \item \code{r_truncated(n, x1, x2)}: Take n draws from distribution g
#' truncated to the open interval \eqn{(x1, x2)}.
#' }
#' 
#' We make use of the \code{fibonacci_heap} structure in the
#' \code{datastructures} package to avoid repeated sorted when
#' \code{small_rects} is used as the method of knot selection.
#' 
#' @references
#' Simon Dirmeier, (2018). datastructures: An R package for organisation and
#' storage of data . Journal of Open Source Software, 3(28), 910,
#' \url{https://doi.org/10.21105/joss.00910}
#' 
#' Winston Chang (2020). R6: Encapsulated Classes with Reference Semantics. R
#' package version 2.5.0. \url{https://CRAN.R-project.org/package=R6}
#' 
#' @export
Stepdown = R6Class("Stepdown",
	lock_objects = FALSE,
	lock_class = FALSE,
	private = list(
		N = NULL,
		tol = NULL,
		log_x_vals = NULL,
		log_h_vals = NULL,
		cum_probs = NULL,
		norm_const = NULL,
		log_p = NULL,
		priority_weight = NULL
	),
	public = list(

	# Note: I couldn't immediately figure out how to get Roxygen working with
	# R6 for functions defined outside of this block using the "set"
	# approach. The added support for R6 in Roxygen seems to be a fairly new
	# development, and support for "set" may be added soon (or is already done
	# and I missed it).

	#' @description Construct the step function
	initialize = function(w, g, tol, N, method, priority_weight = 1/2) {
		stopifnot(class(w) == "weight")
		stopifnot(class(g) == "base")
		private$tol = tol
		private$N = N
		private$priority_weight = priority_weight
		private$setup(w, g, tol, N, method)
	},

	#' @description Get cumulative probabilities after normalizing
	#' the step function to a discrete distribution.
	get_cum_probs = function() {
		private$cum_probs
	},

	#' @description The constant used to normalize the step function
	#' to a discrete distribution.
	get_norm_const = function() {
		private$norm_const
	},

	#' @description The knot values on which the approximation is based,
	#' returned on the log-scale.
	get_log_x_vals = function() {
		private$log_x_vals
	},

	#' @description Values of the approximation evaluated at the points
	#' returned by \code{get_log_x_vals}, given at the log-scale.
	get_log_h_vals = function() {
		private$log_h_vals
	},

	#' @description Evaluate the function \eqn{\rm{Pr}(A_u)} on a given
	#' point u (to be provided on the log-scale). The result is returned
	#' on the log-scale.
	get_log_p = function(log_u) {
		private$log_p(log_u)
	},

	get_priority_weight = function() {
		private$priority_weight
	},

	#' @description Update step function by adding a knot at the point u
	#' (given on the log-scale).
	add = function(log_u)
	{
		# This is not very efficient; we should be able to do better in C++.
		log_h_val = private$log_p(log_u)
		private$log_x_vals = sort(c(private$log_x_vals, log_u))
		private$log_h_vals = sort(c(private$log_h_vals, log_h_val), decreasing = TRUE)
		private$update()
	},

	#' @description Quantiles of the distribution based on the step function.
	q = function(p)
	{
		n = length(p)
		out = numeric(n)
		cum_probs_ext = c(0, private$cum_probs)
		for (i in 1:n) {
			j1 = findInterval(p[i], cum_probs_ext, rightmost.closed = TRUE)
			j2 = j1 + 1
			x1 = exp(private$log_x_vals[j1])
			x2 = exp(private$log_x_vals[j2])
			cp1 = cum_probs_ext[j1]
			cp2 = cum_probs_ext[j2]
			out[i] = x1 + (x2 - x1) * (p - cp1) / (cp2 - cp1)
		}
		return(out)
	},

	#' @description Draw from the distribution based on the step function.
	r = function(n)
	{
		u = runif(n)
		out = numeric(n)
		for (i in 1:n) {
			out[i] = self$q(u[i])
		}
		return(out)
	},

	#' @description Density of the distribution based on the step function.
	d = function(x, log = FALSE, normalize = TRUE)
	{
		n = length(x)
		out = numeric(n)
		for (i in 1:n) {
			# Get the idx such that x_vals[idx] <= x < x_vals[idx+1]
			idx = findInterval(log(x[i]), private$log_x_vals)
			out[i] = private$log_h_vals[idx]
		}
		if (normalize) { out = out - private$norm_const }
		if (log) { return(out) } else { return(exp(out)) }
	},

	#' @description CDF of the distribution based on the step function.
	p = function(x)
	{
		n = length(x)
		out = numeric(n)
		for (i in 1:n) {
			# Get the idx such that x_vals[idx] <= x < x_vals[idx+1]
			j1 = findInterval(log(x[i]), private$log_x_vals)
			j2 = j1 + 1
			x1 = exp(private$log_x_vals[j1])
			x2 = exp(private$log_x_vals[j2])
			cp1 = private$cum_probs[j1]
			cp2 = private$cum_probs[j2]
			out[i] = cp1 + (x[i] - x1) / (x2 - x1) * (cp2 - cp1)
		}
		return(out)
	})
)

# Build a stepdown function
Stepdown$set("private", "setup", function(w, g, tol, N, method)
{
	# The geometric mean
	midpoint = function(log_x, log_y, take_log = TRUE) {
		out = 1/2 * (log_x + log_y)
		ifelse(take_log, out, exp(out))
	}

	# Compute log(y - x) given log(x) and log(y)
	unidist = function(log_x, log_y) { logsub(log_y, log_x) }

	# Compute r*x + (1-r)*y
	# Compute on the log-scale in case we encounter very small magnitude numbers
	log_p = function(log_u) {
		endpoints = w$roots(w$log_c + log_u)
		log(g$pr_interval(endpoints[1], endpoints[2]))
	}

	# Save functions for later use
	private$midpoint = midpoint
	private$log_p = log_p

	# First, make sure A_0 has positive probability according to g.
	# If it doesn't, the rest of the algorithm won't work, so bail out.
	log_prob_max = log_p(-Inf)
	if (is.infinite(log_prob_max)) {
		stop("Could not find any u such that P(X in A_u) > 0")
	}

	# Find L, the largest point where where P(A_L) < prob_max.
	# First find a finite lower bound for log_L
	log_L_lo = -1
	log_L_hi = 0
	log_prob = log_p(log_L_lo)
	while (log_prob < log_prob_max) {
		log_L_lo = 2*log_L_lo
		log_prob = log_p(log_L_lo)
	}

	# Do a bisection search to find log_L between log_L_lo and log_L_hi.
	pred_logL = function(log_u) {
		# Compare prob < prob_max on the log-scale
		log_p(log_u) < log_prob_max
	}
	delta = tol * (log_L_hi - log_L_lo)
	log_L = bisection(log_L_lo, log_L_hi, pred_logL, midpoint, unidist, log_L_lo)

	# Do a bisection search to find U, the smallest point where P(A_U) = 0.	
	# Find the smallest point where P(A_U) > 1e-5 * P(A_L).
	pred_logU = function(log_u) {
		# Compare on the log-scale
		log_p(log_u) < log(1e-5) + log_p(log_L)
	}
	delta = tol * (0 - log_L)
	log_U = bisection(log_L, 0, pred_logU, midpoint, unidist, log_L)

	# Now fill in points between L and U
	if (method == "equal_steps") {
		private$init_equal_steps(log_L, log_U, log_prob_max)
	} else if (method == "small_rects") {
		private$init_small_rects(log_L, log_U, log_prob_max)
	} else {
		msg = sprintf("Unknown method: %s. Currently support %s and %s\n",
          method, "equal_steps", "small_rects")
		stop(msg)
	}

	private$update()
})

Stepdown$set("private", "init_equal_steps", function(log_L, log_U, log_prob_max)
{
	N = private$N
	log_x_vals = log(seq(exp(log_L), exp(log_U), length.out = N+1))
	log_h_vals = numeric(N+1)

	for (i in seq_len(N+1)) {
		log_h_vals[i] = private$log_p(log_x_vals[i])
	}

	# Add pieces for the interval [0, L) where h(u) = 1
	private$log_x_vals = c(-Inf, log_x_vals)
	private$log_h_vals = c(log_prob_max, log_h_vals)

})

Stepdown$set("private", "init_small_rects", function(log_L, log_U, log_prob_max)
{
	N = private$N
	tol = private$tol
	pw = private$priority_weight

	# This queue should be in max-heap order by height
	q = fibonacci_heap("numeric")
	intvl = get_interval(log_L, log_U, private$log_p(log_L), private$log_p(log_U))
	priority = pw * intvl$log_height + (1-pw) * intvl$log_width
	insert(q, -priority, intvl)

	# Preallocate log_x_vals and log_h_vals to the maximum size we will need.
	log_x_vals = rep(NA, N+2)
	log_h_vals = rep(NA, N+2)

	log_x_vals[1] = -Inf
	log_x_vals[2] = log_L
	log_x_vals[3] = log_U

	log_h_vals[1] = log_prob_max
	log_h_vals[2] = private$log_p(log_L)
	log_h_vals[3] = private$log_p(log_U)

	# We already have three (x, h(x)) pairs from above
	iter = 3

	while (size(q) > 0 && iter < N+2) {
		# Get the interval with the highest priority
		int_top = pop(q)[[1]]

		# Break the interval int_top into two pieces: left and right.
		log_x_new = private$midpoint(int_top$log_x, int_top$log_y, take_log = TRUE)
		log_h_new = private$log_p(log_x_new)

		# Add the midpoint to our list of knots
		iter = iter + 1
		log_x_vals[iter] = log_x_new
		log_h_vals[iter] = log_h_new

		# Add the interval which represents [int_top$log_x, log_x_new]
		int_left = get_interval(int_top$log_x, log_x_new, int_top$log_h_x, log_h_new)
		priority = pw * int_left$log_height + (1-pw) * int_left$log_width
		insert(q, -priority, int_left)

		# Add the interval which represents [log_x_new, int_top$log_x]
		int_right = get_interval(log_x_new, int_top$log_y, log_h_new, int_top$log_h_y)
		priority = pw * int_right$log_height + (1-pw) * int_right$log_width
		insert(q, -priority, int_right)
	}

	private$log_x_vals = sort(log_x_vals[seq_len(iter)], decreasing = FALSE)
	private$log_h_vals = sort(log_h_vals[seq_len(iter)], decreasing = TRUE)
})

Stepdown$set("private", "update", function()
{
	N = length(private$log_x_vals) - 2

	# The approximation was computed on the log-scale. Transform it back to the
	# probability sample and integrate over the pieces to find the normalizing
	# constant. We need to be careful here because magnitudes of widths can be
	# extremely tiny floating point numbers.

	if (FALSE) {
		# This is an easier version of the calculation that may be less precise
		widths = diff(exp(private$log_x_vals))
		heights = exp(private$log_h_vals[-(N+2)])
		areas = widths * heights
		private$norm_const = sum(areas)
		private$cum_probs = cumsum(areas) / private$norm_const
	} else {
		# This version tries to be more precise and do more work on the log-scale
		log_areas = numeric(N+1)
		log_cum_areas = numeric(N+1)

		log_areas[1] = private$log_h_vals[1] + private$log_x_vals[2]
		log_cum_areas[1] = log_areas[1]
		for (l in seq(2, N+1)) {
			log_A = private$log_h_vals[l] + private$log_x_vals[l+1]
			log_B = private$log_h_vals[l] + private$log_x_vals[l]
			log_areas[l] = logsub(log_A, log_B)
			log_cum_areas[l] = logadd(log_cum_areas[l-1], log_areas[l])
		}

		log_normconst = log_cum_areas[N+1]
		log_cum_probs = log_cum_areas - log_normconst
		private$norm_const = exp(log_normconst)
		private$cum_probs = exp(log_cum_probs)
	}

	if (any(is.na(private$cum_probs))) {
		warning("NA values found in cum_probs. Approximation may have failed")
	}
})

Stepdown$lock()


# This is used in the heap implementation of init_small_rects.
# It represents an interval [x,y] with function values h(x) and h(y).
get_interval = function(log_x, log_y, log_h_x, log_h_y)
{
	# Consider changing these to have a "log_" prefix
	log_width = logsub(log_y, log_x)
	log_height = logsub(log_h_x, log_h_y)
	ret = list(log_x = log_x, log_y = log_y, log_h_x = log_h_x,
		log_h_y = log_h_y, log_width = log_width, log_height = log_height)
	class(ret) = "interval"
	return(ret)
}

print.interval = function(intvl, log_scale = FALSE)
{
	if (log_scale) {
		printf("log_x: %g\n", intvl$log_x)
		printf("log_y: %g\n", intvl$log_y)
		printf("log_h_x: %g\n", intvl$log_h_x)
		printf("log_h_y: %g\n", intvl$log_h_y)
		printf("width: %g\n", intvl$log_width)
		printf("height: %g\n", intvl$log_height)
	} else {
		printf("x: %g\n", exp(intvl$log_x))
		printf("y: %g\n", exp(intvl$log_y))
		printf("h_x: %g\n", exp(intvl$log_h_x))
		printf("h_y: %g\n", exp(intvl$log_h_y))
		printf("width: %g\n", exp(intvl$log_width))
		printf("height: %g\n", exp(intvl$log_height))
	}
}
