source("gibbs.R")

# These arguments should be set by the driver script
stopifnot(exists("n"))
stopifnot(exists("sigma2_true"))
stopifnot(exists("rho"))
stopifnot(exists("beta_true"))
stopifnot(exists("N_sim"))
stopifnot(exists("mcmc_reps"))
stopifnot(exists("mcmc_burn"))
stopifnot(exists("report_period"))

# ----- Set up data that will be fixed through simulation -----
d = length(beta_true)
X = matrix(NA, n, d)
for (l in 1:d) {
	X[,l] = rnorm(n)
}
mu_true = as.numeric(X %*% beta_true)

prior = get_prior(a_sigma = 2, b_sigma = 10, sigma2_beta = 10^2, d)
control = get_gibbs_control(save_xi = FALSE, R = mcmc_reps,
	burn = mcmc_burn, thin = 1, report_period = report_period,
	direct = get_direct_control(N = 100))

res1 = list()
res2 = list()

for (s in 1:N_sim) {
	logger("*** Simulation %d ***\n", s)

	# Draw new responses and additive noise
	y_true = rlnorm(n, mu_true, sqrt(sigma2_true))
	xi_true = rgeom(n, rho) - rgeom(n, rho)
	y_tilde = y_true + xi_true

	# Fit the noisy observations
	init = get_init(beta = rep(0, d), sigma2 = 1, n = n, d = d)
	fit1_out = my_sampler(y_tilde, X, rho, init, prior, control)

	# Fit observations as if there were no added noise
	init$xi = rep(0,n)
	fit2_out = my_sampler(y_true, X, rho = rep(1,n), init, prior, control,
		fixed = get_fixed(xi = TRUE))

	res1[[s]] = fit1_out
	res2[[s]] = fit2_out
}

save.image("results.Rdata")

