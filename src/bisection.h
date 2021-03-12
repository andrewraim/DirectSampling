#ifndef BISECTION_H
#define BISECTION_H

#include "functionals.h"

double bisection(double x_lo, double x_hi,
	const Predicate& f, const Functional2& mid,
	const Functional2& dist, double tol);

#endif
