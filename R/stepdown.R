# This is used in the heap implementation of init_small_rects.
# It represents an interval [x,y] with function values h(x) and h(y).
# We save width = y-x and height = log[ h(x) / h(y) ].
get_interval = function(log_x, log_y, log_h_x, log_h_y)
{
	# Consider changing these to have a "log_" prefix
	width = exp(log_y) - exp(log_x)
	height = exp(log_h_x) - exp(log_h_y)
	priority = log(height) + log(width)
	ret = list(log_x = log_x, log_y = log_y, log_h_x = log_h_x,
		log_h_y = log_h_y, width = width, height = height,
		priority = priority)
	class(ret) = "interval"
	return(ret)
}

print.interval = function(intvl)
{
	printf("log_x: %g\n", intvl$log_x)
	printf("log_y: %g\n", intvl$log_y)
	printf("log_h_x: %g\n", intvl$log_h_x)
	printf("log_h_y: %g\n", intvl$log_h_y)
	printf("width: %g\n", intvl$width)
	printf("height: %g\n", intvl$height)
	printf("priority: %g\n", intvl$priority)
}

# This R6 class represents a step function approximation to the non-increasing
# Pr(A_u) for u in [0,1].
#' @export
Stepdown = R6Class("Stepdown",
	lock_objects = FALSE,
	lock_class = FALSE,
	private = list(
		log_x_vals = NULL,
		log_h_vals = NULL,
		cum_probs = NULL,
		norm_const = NULL,
		log_p = NULL
	),
	public = list(
		initialize = function(w, g, tol, N, method) {
			stopifnot(class(w) == "weight")
			stopifnot(class(g) == "base")
			private$tol = tol
			private$N = N
			private$setup(w, g, tol, N, method)
		},
		get_cum_probs = function() {
			private$cum_probs
		},
		get_norm_const = function() {
			private$norm_const
		},
		get_log_x_vals = function() {
			private$log_x_vals
		},
		get_log_h_vals = function() {
			private$log_h_vals
		},
		get_log_p = function(log_u) {
			private$log_p(log_u)
		}
	)
)

# Build a stepdown function
Stepdown$set("private", "setup", function(w, g, tol, N, method)
{
	midpoint = function(x,y) { (x+y)/2 }
	unidist = function(x,y) { y-x }
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
	log_L = bisection(log_L_lo, log_L_hi, pred_logL, midpoint, unidist, delta)

	# Do a bisection search to find U, the largest point where P(A_U) > 0.	
	pred_logU = function(log_u) {
		# Compare on the log-scale
		is.infinite(log_p(log_u))
	}
	delta = tol * (0 - log_L)
	log_U = bisection(log_L, 0, pred_logU, midpoint, unidist, delta)

	# Now fill in points between L and U
	if (method == "equal_steps") {
		private$init_equal_steps(log_L, log_U, log_prob_max, N)
	} else if (method == "small_rects") {
		private$init_small_rects(log_L, log_U, log_prob_max, N)
	} else {
		msg = sprintf("Unknown method: %s. Currently support %s and %s\n",
          method, "equal_steps", "small_rects")
		stop(msg)
	}

	private$update()
})

Stepdown$set("private", "init_equal_steps", function(log_L, log_U,
	log_prob_max, N)
{
	log_x_vals = log(seq(exp(log_L), exp(log_U), length.out = N+1))
	log_h_vals = numeric(N+1)

	for (i in seq_len(N+1)) {
		log_h_vals[i] = private$log_p(log_x_vals[i])
	}

	# Add pieces for the interval [0, L) where h(u) = 1
	private$log_x_vals = c(-Inf, log_x_vals)
	private$log_h_vals = c(log_prob_max, log_h_vals)

})

Stepdown$set("private", "init_small_rects", function(log_L, log_U,
	log_prob_max, N)
{
	# This queue should be in max-heap order by height
	q = fibonacci_heap("numeric")
	intvl = get_interval(log_L, log_U, private$log_p(log_L), private$log_p(log_U))
	insert(q, -intvl$priority, intvl)

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
		# Get the interval with the largest height
		int_top = pop(q)[[1]]

		# Break the interval int_top into two pieces: left and right.
		log_x_new = log(private$midpoint(exp(int_top$log_x), exp(int_top$log_y)))
		log_h_new = private$log_p(log_x_new)

		# Add the midpoint to our list of knots
		iter = iter + 1
		log_x_vals[iter] = log_x_new
		log_h_vals[iter] = log_h_new

		# Add the interval which represents [int_top$log_x, log_x_new]
		int_left = get_interval(int_top$log_x, log_x_new, int_top$log_h_x, log_h_new)
		insert(q, -int_left$priority, int_left)

		# Add the interval which represents [log_x_new, int_top$log_x]
		int_right = get_interval(log_x_new, int_top$log_y, log_h_new, int_top$log_h_y)
		insert(q, -int_right$priority, int_right)
	}

	private$log_x_vals = sort(log_x_vals[seq_len(iter)], decreasing = FALSE)
	private$log_h_vals = sort(log_h_vals[seq_len(iter)], decreasing = TRUE)
})

Stepdown$set("public", "q", function(p)
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
})

Stepdown$set("public", "r", function(n)
{
	u = runif(n)
	out = numeric(n)
	for (i in 1:n) {
		out[i] = self$q(u[i])
	}
	return(out)
})

Stepdown$set("public", "d", function(x, log = FALSE, normalize = TRUE)
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
})

Stepdown$set("public", "p", function(x)
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

Stepdown$set("public", "add", function(log_u)
{
	# Update step function with new log_u and log_p values.
	# This is not very efficient; we will do better in C++.
	log_h_val = private$log_p(log_u)
	private$log_x_vals = sort(c(private$log_x_vals, log_u))
	private$log_h_vals = sort(c(private$log_h_vals, log_h_val), decreasing = TRUE)
	private$update()
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
	} else {
		# This version tries to be more precise and do more work on the log-scale
		idx1 = 1:(N+1)
		idx2 = 2:(N+2)
		areas = numeric(N+1)
		areas[1] = exp(private$log_h_vals[1] + private$log_x_vals[2])
		areas[idx1] = exp(private$log_h_vals[idx1] + private$log_x_vals[idx2]) -
			exp(private$log_h_vals[idx1] + private$log_x_vals[idx1])
	}

	private$norm_const = sum(areas)
	private$cum_probs = cumsum(areas) / private$norm_const

	if (any(is.na(private$cum_probs))) {
		warning("NA values found in cum_probs. Approximation may have failed")
	}
})

Stepdown$lock()

