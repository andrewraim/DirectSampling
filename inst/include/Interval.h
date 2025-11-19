#ifndef DIRECT_SAMPLING_INTERVAL_H
#define DIRECT_SAMPLING_INTERVAL_H

#include <Rcpp.h>
#include "log-sum-exp.h"

namespace DirectSampling {

/*
* This class is used in the heap implementation of Stepdown::init_small_rects.
* It represents an interval [x,y] with function values h(x) and h(y).
* We save width = y-x and height = h(x) - h(y).
*/
class Interval
{
public:
	Interval(double log_x, double log_y, double log_h_x, double log_h_y, double pw = 0.5)
		: _log_x(log_x), _log_y(log_y), _log_h_x(log_h_x), _log_h_y(log_h_y), _pw(pw)
	{
		if (_pw < 0 || _pw > 1) {
			Rcpp::stop("pw must be in [0,1]");
		}
	}
	Interval(const Interval& rhs) {
		*this = rhs;
	}
	double log_width() const {
		return log_sub2_exp(_log_y, _log_x);
	}
	double log_height() const {
		return log_sub2_exp(_log_h_x, _log_h_y);
	}
	double priority() const {
		double log_width = log_sub2_exp(_log_y, _log_x);
		double log_height = log_sub2_exp(_log_h_x, _log_h_y);
		return _pw * log_height + (1-_pw) * log_width;
	}
	const Interval& operator=(const Interval& rhs) {
		_log_x = rhs._log_x;
		_log_y = rhs._log_y;
		_log_h_x = rhs._log_h_x;
		_log_h_y = rhs._log_h_y;
		_pw = rhs._pw;
		return *this;
	}
	bool operator<(const Interval& rhs) const {
		return priority() < rhs.priority();
	}
	bool operator>(const Interval& rhs) const {
		return priority() > rhs.priority();
	}
	double get_log_x() const {
		return _log_x;
	}
	double get_log_y() const {
		return _log_y;
	}
	double get_log_h_x() const {
		return _log_h_x;
	}
	double get_log_h_y() const {
		return _log_h_y;
	}
	void print(bool log_scale = false) const {
		if (log_scale) {
			Rprintf("log_x: %g\n", _log_x);
			Rprintf("log_y: %g\n", _log_y);
			Rprintf("log_h_x: %g\n", _log_h_x);
			Rprintf("log_h_y: %g\n", _log_h_y);
			Rprintf("width: %g\n", log_width());
			Rprintf("height: %g\n", log_height());
			Rprintf("priority: %g\n", priority());
		} else {
			Rprintf("x: %g\n", exp(_log_x));
			Rprintf("y: %g\n", exp(_log_y));
			Rprintf("h_x: %g\n", exp(_log_h_x));
			Rprintf("h_y: %g\n", exp(_log_h_y));
			Rprintf("width: %g\n", exp(log_width()));
			Rprintf("height: %g\n", exp(log_height()));
			Rprintf("priority: %g\n", priority());
		}
	}
private:
	double _log_x;
	double _log_y;
	double _log_h_x;
	double _log_h_y;
	double _pw;
};

}

#endif
