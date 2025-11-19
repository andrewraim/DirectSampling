#ifndef DIRECT_SAMPLING_REJECTION_ARGS_H
#define DIRECT_SAMPLING_REJECTION_ARGS_H

#include "fntl.h"

namespace DirectSampling {

struct rejection_args
{
	unsigned int max_rejects = std::numeric_limits<unsigned int>::max();
	unsigned int report = std::numeric_limits<unsigned int>::max();
	bool adapt = true;	
	fntl::error_action action = fntl::error_action::STOP;

	rejection_args() { };
};

}

#endif

