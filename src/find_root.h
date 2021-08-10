#ifndef FIND_ROOT_H
#define FIND_ROOT_H

#include <Rcpp.h>
#include "bisection.h"

double find_root(double x_lo, double x_hi, const Functional1& f, double tol);

#endif
