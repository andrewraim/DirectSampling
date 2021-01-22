# x_lo: lower bound for solution, should be finite
# x_hi: upper bound for solution, should be finite
# pred(x): predicate function
# mid(x,y): a midpoint function which returns a point between x <= y
# dist(x,y): a a distance functions between x <= y
# tol: Stop when dist(x_lo,x_hi) < tol
# Watch out for numerical overflow. eventually mid(x,y) will no longer be
# strictly between x and y
#' @export
bisection = function(x_lo, x_hi, pred, mid, dist, tol)
{
	x = mid(x_lo, x_hi)

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi) {
		ind = pred(x)
		# printf("ind = %d, x_lo = %24.16f, x_hi = %24.16f\n", ind, x_lo, x_hi)
		x_lo = ind*x_lo + (1-ind)*x
		x_hi = ind*x + (1-ind)*x_hi
		x = mid(x_lo, x_hi)
	}

	if (dist(x_lo, x_hi) > tol && (x <= x_lo || x >= x_hi)) {
		stop("Numerical overflow in bisection. tol may be too small")
	}

	return(x)
}

