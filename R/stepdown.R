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
#' @param log_x A vector of points to evaluate, provided on the log-scale.
#' @param log If \code{TRUE}, return value on the log-scale.
#' @param normalize If \code{TRUE}, normalize the result to be a density value.
#' @param priority_weight Experimental: An weight between 0 and 1. When closer
#' to 1, more priority is given in \code{small_rects} to rectangle height. When
#' closer to 0, more priority is given to rectangle width.
#' @param midpoint_type Type of midpoint function to use. Currently only
#' \code{geometric} or \code{arithmetic} are supported. \code{geometric} is
#' preferred when the support of \eqn{p(u)} is focused very close to zero.
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
		knot_order = NULL,
		cum_probs = NULL,
		norm_const = NULL,
		log_p = NULL,
		priority_weight = NULL,
		midpoint_type = NULL
	),
	public = list(

	# Note: I couldn't immediately figure out how to get Roxygen working with
	# R6 for functions defined outside of this block using the "set"
	# approach. The added support for R6 in Roxygen seems to be a fairly new
	# development, and support for "set" may be added soon (or is already done
	# and I missed it).

	#' @description Construct the step function
	initialize = function(w, g, tol, N, method, priority_weight = 1/2,
		midpoint_type = "geometric") {
		stopifnot(class(w) == "weight")
		stopifnot(class(g) == "base")
		private$tol = tol
		private$N = N
		private$priority_weight = priority_weight
		private$midpoint_type = midpoint_type
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

	#' @description Rectangles used to bound the approximation: x and h
	#' coordinates along with the area. Values are returned on the log-scale or
	#' the original scale according to the \code{log} argument.
	get_rects = function(log = FALSE) {
		k = private$N + 2
		log_x1 = private$log_x_vals[-k]
		log_x2 = private$log_x_vals[-1]
		log_h1 = private$log_h_vals[-k]
		log_h2 = private$log_h_vals[-1]
		log_area = numeric(k-1)
		for (i in seq_len(k-1)) {
			log_area[i] = logsub(log_x2[i], log_x1[i]) +  logsub(log_h1[i], log_h2[i])
		}

		if (log) {
			rects = data.frame(log_x1, log_x2, log_h1, log_h2, log_area)
		} else {
			rects = data.frame(
				x1 = exp(log_x1),
				x2 = exp(log_x2),
				h1 = exp(log_h1),
				h2 = exp(log_h2),
				area = exp(log_area)
			)
		}
		return(rects)
	},

	#' @description Returns the order in which knot points were added to the
	#' step function. Indices correspond to the points returned by
	#' \code{get_log_x_vals}.
	get_knot_order = function() {
		private$knot_order
	},

	#' @description Evaluate the function \eqn{\rm{Pr}(A_u)} on a given
	#' point u (to be provided on the log-scale). The result is returned
	#' on the log-scale.
	get_log_p = function(log_u) {
		private$log_p(log_u)
	},

	#' @description Returns the value of \code{priority_weight}.
	get_priority_weight = function() {
		private$priority_weight
	},

	#' @description Update step function by adding a knot at the point u
	#' (given on the log-scale).
	add = function(log_u)
	{
		log_h = private$log_p(log_u)

		# Add new x and h value, and keep track of the order in which it was added
		# This is not very efficient; we could potentially to do better than resorting.
		log_x_vals = c(private$log_x_vals, log_u)
		log_h_vals = c(private$log_h_vals, log_h)
		knot_order = c(private$knot_order, private$N + 3)
		idx = order(log_x_vals)
		private$log_x_vals = log_x_vals[idx]
		private$log_h_vals = log_h_vals[idx]
		private$knot_order = knot_order[idx]

		# Update cumulative probabilities
		private$N = private$N + 1
		private$update()
	},

	#' @description Quantiles of the distribution based on the step function.
	q = function(p, log = FALSE)
	{
		n = length(p)
		out = numeric(n)
		cum_probs_ext = c(0, private$cum_probs)
		for (i in 1:n) {
			j1 = findInterval(p[i], cum_probs_ext, rightmost.closed = TRUE)
			j2 = j1 + 1
			cp1 = cum_probs_ext[j1]
			cp2 = cum_probs_ext[j2]
			log_x1 = private$log_x_vals[j1]
			log_x2 = private$log_x_vals[j2]

			# Do the following computation, but be careful to keep x values on
			# the log-scale, since they may be extremely small.
			# out[i] = x1 + (x2 - x1) * (p[i] - cp1) / (cp2 - cp1)
			log_ratio = log(p[i] - cp1) - log(cp2 - cp1)
			out[i] = logadd(logsub(log_x2, log_x1) + log_ratio, log_x1)
		}
		if (any(is.na(out))) {
			stop("NA or NaN values were produced in Stepdown::q")
		}

		if (log) { return(out) } else { return(exp(out)) }
	},

	#' @description Draw from the distribution based on the step function.
	r = function(n, log = TRUE)
	{
		u = runif(n)
		out = numeric(n)
		for (i in 1:n) {
			out[i] = self$q(u[i], log = log)
		}
		return(out)
	},

	#' @description Density of the distribution based on the step function.
	d = function(log_x, log = FALSE, normalize = TRUE)
	{
		n = length(log_x)
		out = rep(-Inf, n)
		for (i in 1:n) {
			# Get the idx such that x_vals[idx] <= x < x_vals[idx+1]
			idx = findInterval(log_x[i], private$log_x_vals)
			if (idx < private$N + 2) {
				out[i] = private$log_h_vals[idx]
			}
		}
		if (normalize) { out = out - private$norm_const }
		if (log) { return(out) } else { return(exp(out)) }
	},

	#' @description CDF of the distribution based on the step function.
	p = function(log_x)
	{
		n = length(log_x)
		out = numeric(n)
		for (i in 1:n) {
			# Get the idx such that x_vals[idx] <= x < x_vals[idx+1]
			j1 = findInterval(log_x[i], private$log_x_vals)
			j2 = j1 + 1
			cp1 = private$cum_probs[j1]
			cp2 = private$cum_probs[j2]
			log_x1 = private$log_x_vals[j1]
			log_x2 = private$log_x_vals[j2]

			# Do the following computation, but be careful to keep x values on
			# the log-scale, since they may be extremely small.
			# out[i] = cp1 + (x[i] - x1) / (x2 - x1) * (cp2 - cp1)
			log_ratio = logsub(log_x[i], log_x1) - logsub(log_x2, log_x1) + log(cp2 - cp1)
			out[i] = exp(logadd(log(cp1), log_ratio))
		}
		if (any(is.na(out))) {
			stop("NA or NaN values were produced in Stepdown::q")
		}
		return(out)
	})
)

# Build a stepdown function
Stepdown$set("private", "setup", function(w, g, tol, N, method)
{
	if (private$midpoint_type == "geometric") {
		# The geometric mean
		midpoint = function(log_x, log_y, take_log = TRUE) {
			out = 1/2 * (log_x + log_y)
			ifelse(take_log, out, exp(out))
		}
	} else if (private$midpoint_type == "arithmetic") {
		# The arithmetic mean, computed on the log-scale
		midpoint = function(log_x, log_y, take_log = TRUE) {
			out = log(1/2) + log_y + log1p(exp(log_x - log_y))
			ifelse(take_log, out, exp(out))
		}
	} else {
		stop("Unrecognized midpoint type")
	}

	# Compute log(y - x) given log(x) and log(y)
	unidist = function(log_x, log_y) { logsub(log_y, log_x) }

	# Compute on the log-scale in case we encounter very small magnitude numbers
	log_p = function(log_u) {
		endpoints = w$roots(w$log_c + log_u)
		g$pr_interval(endpoints[1], endpoints[2], log = TRUE)
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
	log_delta1 = min(log_L_lo, log(tol))
	log_L = bisection(log_L_lo, log_L_hi, pred_logL, midpoint, unidist, log_delta1)

	# Do a bisection search to find U, the smallest point where P(A_U) = 0.	
	# Find the smallest point where P(A_U) > small_factor * P(A_L).
	pred_logU = function(log_u) {
		# Compare on the log-scale
		# TBD: small_factor = 1e-10 is a magic number for now. If it's too small,
		# we'll get areas where nothing much is happening in the support. But if
		# it's too large, we might miss a little bit of useful support.
		log_p(log_u) < log(1e-10) + log_p(log_L)
	}
	log_delta2 = min(log_L, log(tol))
	log_U = bisection(log_L, 0, pred_logU, midpoint, unidist, log_delta2)

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
	private$knot_order = seq_len(N+2)

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
		log_x_new = private$midpoint(int_top$log_x, int_top$log_y)
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

	idx = order(log_x_vals)
	private$log_x_vals = log_x_vals[idx]
	private$log_h_vals = log_h_vals[idx]
	private$knot_order = idx
})

Stepdown$set("private", "update", function()
{
	N = private$N

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
			log_areas[l] = private$log_h_vals[l] + logsub(private$log_x_vals[l+1], private$log_x_vals[l])
			arg1 = max(log_cum_areas[l-1], log_areas[l])
			arg2 = min(log_cum_areas[l-1], log_areas[l])
			log_cum_areas[l] = logadd(arg1, arg2)
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
	if (log_h_x < log_h_y) {
		log_height = 0
	} else {
		log_height = logsub(log_h_x, log_h_y)
	}
	ret = list(log_x = log_x, log_y = log_y, log_h_x = log_h_x,
		log_h_y = log_h_y, log_width = log_width, log_height = log_height)
	class(ret) = "interval"
	return(ret)
}

#' @export
print.interval = function(x, log_scale = FALSE, ...)
{
	if (log_scale) {
		printf("log_x: %g\n", x$log_x)
		printf("log_y: %g\n", x$log_y)
		printf("log_h_x: %g\n", x$log_h_x)
		printf("log_h_y: %g\n", x$log_h_y)
		printf("width: %g\n", x$log_width)
		printf("height: %g\n", x$log_height)
	} else {
		printf("x: %g\n", exp(x$log_x))
		printf("y: %g\n", exp(x$log_y))
		printf("h_x: %g\n", exp(x$log_h_x))
		printf("h_y: %g\n", exp(x$log_h_y))
		printf("width: %g\n", exp(x$log_width))
		printf("height: %g\n", exp(x$log_height))
	}
}

