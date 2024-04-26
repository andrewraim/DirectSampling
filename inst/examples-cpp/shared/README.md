---
title: >
  Replication Materials for Manuscript "Direct Sampling with a Step Function":
  Rcpp folder
author: Andrew Raim
---

C++ code used in the examples. Requires Rcpp to compile and run.

- direct_sampler.{h,cpp}: direct sampler with step function but without
  accept-reject algorithm.
- direct_sampler_ar.{h,cpp}: direct sampler with step function and also with
  accept-reject.
- Stepdown.{h,cpp}: a class that represents the step function. This is the most
  involved part of the code.
- bisection.{h,cpp}: the bisection algorithm described in the manuscript.
- find_interval.{h,cpp}: code to find an interval which contains the specified
  value.
- find_root.{h,cpp}: find the root of a given functional.
- functionals.h: defines a few functional types that are used in the code.
- util.{h,cpp}: some utilities used in the code.
- BaseDistribution.h, WeightFunction.h: abstract base classes for base
  distributions and weight functions, which are notated `g` and `w` in the
  manuscript, respectively.
- The following are specific BaseDistributions and WeightFunctions used in
  examples.
	+ BetaBaseDistribution.h
	+ UniformBaseDistribution.h
	+ CARWeightFunction.h
	+ TDistDFWeightFunction.h

