---
title: >
  Companion codes for the manuscript "Direct Sampling in Bayesian Regression
  Models with Additive Disclosure Avoidance Noise"
author: Andrew Raim
output: 
  html_document:
    number_sections: true
---

# Requirements

This code requires R with the DirectSampling package. Additionally, the
following R packages are used.

- Rcpp
- R6
- ggplot2
- datastructures
- Matrix
- dplyr
- gridExtra

# Contents
To run any of the scripts below, ensure that your current working directory
(e.g. `getwd()`) matches the directory of the script.

## `01-direct-sampler/`
Code for Sections 3 and 4 of the manuscript.

1. `functions.R` contains supporting R functions used by other scripts in this
   folder. It is not intended to be executed directly.

2. The following scripts produce plots shown in Fig 1.

   a. `step1.R`: produce Figs 1c and 1d.
   b. `step2.R`: produce Figs 1a and 1b.

3. The following scripts demonstrate specific examples of direct samplers with
   the specified weight functions and base distributions.

   a. `lognormal-laplace.R`: Lognormal weight and Laplace base.
   b. `lognormal-dgeom.R`: Lognormal weight and DGeom base.
   c. `normal-laplace.R`: Normal weight and Laplace base.
   d. `normal-dgeom.R`:  Normal weight and DGeom base.

## `02-regression/`
Code for Section 5.1 of the manuscript.

1. The file `gibbs.R` implements the Gibbs sampler described by Algorithm 4.
   Algorithm 5 can also be obtained by setting the arguments appropriately.

2. To run the Gibbs sampler on the individual datasets used to produce Figs 4
   and 5, use the script `run-one.R`.

3. Several scripts are provided to generate the simulation directory structure
   carry a run of the simulation, and post-process the results.

   a. `setup.R`: Generate directory structure for the simulation.
   b. `sim.R`: Run one setting of the simulation whose results are shown in
      Fig 6.
   c. `analyze.R`: analyze results from simulation. It assumes that results are
      laid out in the directory structure specified by `setup.R`.

## `03-regression-xnoise/`
Code for Section 5.2 of the manuscript.

1. The file `gibbs.R` implements the Gibbs sampler described by Algorithm 3.
   Algorithm 5 can also be obtained by setting the arguments appropriately.

2. To run the Gibbs sampler on the individual datasets used to produce Figs 7
   and 8, use the script `run-one.R`.

3. Several scripts are provided to generate the simulation directory structure
   carry a run of the simulation, and post-process the results.

   a. `setup.R`: Generate directory structure for the simulation.
   b. `sim.R`: Run one setting of the simulation whose results are shown in
      Fig 9.
   c. `analyze.R`: analyze results from simulation. It assumes that results are
      laid out in the directory structure specified by `setup.R`.


## `shared/`
Code used in several places throughout the other folders.

1. `plots.R`: Code for plotting.
2. `util.R`: Utility functions.

