#ifndef DIRECT_SAMPLING_FUNCTIONALS_H
#define DIRECT_SAMPLING_FUNCTIONALS_H

namespace DirectSampling {

/*
* Typedefs
*/
typedef std::function<bool(double)> predicate_t;
typedef std::function<double(double a, double b)> midpoint_t;
typedef std::function<double(double a, double b)> distance_t;

}

#endif
