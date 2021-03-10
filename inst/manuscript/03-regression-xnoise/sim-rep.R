source("gibbs.R")

# These arguments should be set by the driver script
stopifnot(exists("n"))
stopifnot(exists("sigma2_true"))
stopifnot(exists("rho_y"))
stopifnot(exists("lambda_x"))
stopifnot(exists("beta_true"))
stopifnot(exists("N_sim"))
stopifnot(exists("mcmc_reps"))
stopifnot(exists("mcmc_burn"))
stopifnot(exists("N_direct_y"))
stopifnot(exists("N_direct_x"))
stopifnot(exists("report_period"))

# ----- Set up data that will be fixed through simulation -----
d = length(beta_true)
X_true = matrix(NA, n, d)
for (l in 1:d) {
	X_true[,l] = rnorm(n)
}
mu_true = as.numeric(X_true %*% beta_true)

# Setup for both samplers
prior = get_prior(a_sigma = 2, b_sigma = 10, sigma2_beta = 10^2, d)
control = get_gibbs_control(save_xi_y = FALSE, save_xi_x = FALSE,
	R = mcmc_reps, burn = mcmc_burn, thin = 1, report_period = report_period,
	direct_y = get_direct_control(N = N_direct_y),
	direct_x = get_direct_control(N = N_direct_x))

res1 = list()
res2 = list()

for (s in 1:N_sim) {
	logger("*** Rep %d ***\n", s)

	# Draw new responses and additive noise
	y_true = rlnorm(n, mu_true, sqrt(sigma2_true))
	xi_y_true = rgeom(n, rho_y) - rgeom(n, rho_y)
	y_tilde = y_true + xi_y_true

	# Draw noise for covariates
	xi_x_true = matrix(0, n, d)
	for (i in 1:n) {
		for (j in 1:d1) {
			xi_x_true[i,j] = r_laplace(1, 0, lambda_x[i,j])
		}
	}
	X_tilde = X_true + xi_x_true

	logger("Fitting noisy obs\n")
	init = get_init(n = n, d = d)
	fit1_out = my_sampler(y_tilde, X_tilde, rho_y, lambda_x, init, prior, control)

	logger("Fitting noiseless obs\n")
	init = get_init(xi_y = rep(0,n), xi_x = matrix(0,n,d), n = n, d = d)
	fixed = get_fixed(xi_y = TRUE, xi_x = TRUE)
	fit2_out = my_sampler(y_true, X_true, rho_y = rep(1,n),
		lambda_x = matrix(0,n,d1), init, prior, control, fixed)

	res1[[s]] = fit1_out
	res2[[s]] = fit2_out
}
