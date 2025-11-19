#ifndef DIRECT_SAMPLING_STEPDOWN_ARGS_H
#define DIRECT_SAMPLING_STEPDOWN_ARGS_H

#include "fntl.h"

namespace DirectSampling {

struct stepdown_args
{
	double tol = 1e-8;
	unsigned int N = 30;
	const std::string method = "small_rects";
	double priority_weight = 0.5;

	stepdown_args() { };
};

}

#endif

