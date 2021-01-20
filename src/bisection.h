#ifndef BISECTION_H
#define BISECTION_H

#include <Rcpp.h>
#include "functionals.h"

double bisection(double x_lo, double x_hi,
	const Predicate& pred,
	const Functional2& mid,
	const Functional2& dist, double tol);

#endif
