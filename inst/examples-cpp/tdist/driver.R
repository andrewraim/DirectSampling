library(DirectSampling)
library(tidyverse)
library(splines)

source("gibbs.R")

mcmc_reps = 10000
mcmc_burn = 5000
mcmc_thin = 1
mcmc_report = 1000

# ----- Load data -----
rmbf = read.csv("rmbf.csv")
x = rmbf$x
y = rmbf$y
n = length(y)

# ----- Run Gibbs sampler with direct sampler -----
set.seed(123)

# Design matrix from basis function
X = model.matrix(~ bs(x))

prior = get_prior(nu_min = 0.01, nu_max = 200, a_sigma = 1, b_sigma = 1, sig2beta = 100)
direct_control = get_direct_control(N = 30, tol = 1e-10,
	fill_method = "small_rects", max_rejections = 1000)
control = get_gibbs_control(save_s = FALSE, R = mcmc_reps, burn = mcmc_burn,
	thin = mcmc_thin, report_period = mcmc_report, direct = direct_control)

init = get_init(n = n, d = ncol(X))
gibbs_out = gibbs_sampler(y, X, init, prior, control, fixed = get_fixed())
print(gibbs_out)
