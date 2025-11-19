#ifndef DIRECT_SAMPLING_LOG_SUM_EXP_H
#define DIRECT_SAMPLING_LOG_SUM_EXP_H

#include <Rcpp.h>

namespace DirectSampling {

/*
*  The following macros help us to use the properties:
*  log(x + y) = log(x) + log(1 + y/x)
*             = log(x) + log1p(exp(log(y) - log(x))),
*  log(x - y) = log(x) + log(1 - y/x)
*             = log(x) + log1p(-exp(log(y) - log(x))),
*  with log(x) and log(y) given as inputs, i.e. on the log-scale. When x and y
*  are of very different magnitudes, this is more stable when x is taken to be
*  the larger of the inputs. The most extreme case is when one of the inputs
*  might be -inf; in this case that input should be the second one.
*
*  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
*/

inline double log_add2_exp(double x, double y)
{
	return x + log1p(exp(y - x));
}

inline double log_sub2_exp(double x, double y)
{
	return x + log1p(-exp(y - x));
}

}

#endif
