#include "DirectSampling.h"
#include "util.h"

Rcpp::StringVector timestamp_rcpp()
{
	return DirectSampling::timestamp();
}
