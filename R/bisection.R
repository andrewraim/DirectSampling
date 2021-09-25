#' Bisection Algorithm
#' 
#' Given an increasing step function \eqn{f(x)}, locate the point where the
#' step occurs.
#' 
#' @param x_lo A finite lower bound for \eqn{x} which should occur before the step.
#' @param x_hi A finite upper bound for \eqn{x} which should occur after the step.
#' @param f The function \eqn{f(x)}.
#' @param mid A midpoint function which operates on two points \code{(x,y)} where \code{x <= y}.
#' @param dist A distance function which operates on two points \code{(x,y)} where \code{x <= y}.
#' @param tol Stop when \code{dist(x_lo, x_hi) < tol}.
#'
#' @export
bisection = function(x_lo, x_hi, f, mid, dist, tol)
{
	x = mid(x_lo, x_hi)

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		ind = f(x)
		x_lo = ind*x_lo + (1-ind)*x
		x_hi = ind*x + (1-ind)*x_hi
		x = mid(x_lo, x_hi)
		# printf("mid(%g, %g) = %g\n", x_lo, x_hi, x)
	}

	# Watch out for numerical overflow. eventually mid(x,y) will no longer be
	# strictly between x and y.
	# if (dist(x_lo, x_hi) > tol && (x <= x_lo || x >= x_hi)) {
	# 	stop("Numerical overflow in bisection. tol may be too small")
	# }

	return(x)
}
