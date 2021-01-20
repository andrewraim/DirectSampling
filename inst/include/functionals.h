#ifndef FUNCTIONALS_H
#define FUNCTIONALS_H

/*
* Some simple functional types are needed for the bisection algorithm.
* Functional code can be facilitated with the <functional> library in STL
* but we want this code to work with older compilers.
*/

class Predicate
{
public:
	virtual bool operator()(double x) const = 0;
};

class Functional2
{
public:
	virtual double operator()(double x, double y) const = 0;
};

/*
* We also define a few specific functionals that are reused in several places
* - The length of an interval
* - The midpoint of an interval
*/
class IntervalLength : public Functional2
{
public:
	virtual double operator()(double x, double y) const
	{
		return y - x;
	}
};

class Midpoint : public Functional2
{
public:
	virtual double operator()(double x, double y) const
	{
		return (y + x) / 2;
	}
};

#endif

