library(ggplot2)
source("gibbs.R")
source("../shared/plots.R")

set.seed(1234)

n = 200
d1 = 1
rho_levels = c(0.01, 0.1)
lambda_levels = c(0.05, 0.20)
sigma2_true = 1^2
beta_true = c(5, -1)

x_seq = seq(-400,400)
p_seq = d_dgeom(x_seq, rho_levels[1], take_log = FALSE)
dat_plot = data.frame(x = x_seq, p = p_seq)
g = ggplot(dat_plot, aes(x = x, y = p)) +
	geom_point() +
	theme_bw()
plot(g)

x_seq = seq(-50,50)
p_seq = d_dgeom(x_seq, rho_levels[2], take_log = FALSE)
dat_plot = data.frame(x = x_seq, p = p_seq)
g = ggplot(dat_plot, aes(x = x, y = p)) +
	geom_point() +
	theme_bw()
plot(g)

# ----- Generate data from regression model -----
d = length(beta_true)
X_true = matrix(NA, n, d)
for (l in 1:d) {
	X_true[,l] = rnorm(n)
}

mu_true = X_true %*% beta_true
y_true = rlnorm(n, mu_true, sqrt(sigma2_true))

res_noisy = list()
res_noisefree = list()
idx = 0

# ----- Fit one simulated dataset at each setting -----
for (idx_rho in seq_along(rho_levels)) {
for (idx_lambda in seq_along(lambda_levels)) {
	rho_y = rep(rho_levels[idx_rho], n)
	lambda_x = matrix(lambda_levels[idx_lambda], n, d1)

	#  Generate data from the assumed model
	xi_y_true = rgeom(n, rho_y) - rgeom(n, rho_y)
	y_tilde = y_true + xi_y_true

	xi_x_true = matrix(0, n, d)
	for (i in 1:n) {
		for (j in 1:d1) {
			xi_x_true[i,j] = r_laplace(1, 0, lambda_x[i,j])
		}
	}
	X_tilde = X_true + xi_x_true

	# Setup for both samplers
	prior = get_prior(a_sigma = 2, b_sigma = 10, sigma2_beta = 10^2, d)
	control = get_gibbs_control(save_xi_y = TRUE, save_xi_x = TRUE,
		R = 2000, burn = 1000, report_period = 100,
		direct_y = get_direct_control(N = 100),
		direct_x = get_direct_control(N = 100))

	logger("Fitting noisy obs\n")
	init = get_init(n = n, d = d)
	sampler_out = my_sampler(y_tilde, X_tilde, rho_y, lambda_x, init, prior, control)

	logger("Fitting noiseless obs\n")
	init = get_init(xi_y = rep(0,n), xi_x = matrix(0,n,d), n = n, d = d)
	fixed = get_fixed(xi_y = TRUE, xi_x = TRUE)
	noisefree_out = my_sampler(y_true, X_true, rho_y = rep(1,n),
		lambda_x = matrix(0,n,d1), init, prior, control, fixed)

	idx = idx + 1
	res_noisy[[idx]] = sampler_out
	res_noisefree[[idx]] = noisefree_out
}
}

# ----- Plot results -----

# Common x- and y-limits for all the beta plots
x_lo = min(
	unlist(Map(function(x) { min(x$beta_hist[,1]) }, res_noisy)),
	unlist(Map(function(x) { min(x$beta_hist[,1]) }, res_noisefree))
)
x_hi = max(
	unlist(Map(function(x) { max(x$beta_hist[,1]) }, res_noisy)),
	unlist(Map(function(x) { max(x$beta_hist[,1]) }, res_noisefree))
)
x_ran = c(x_lo, x_hi)

y_lo = min(
	unlist(Map(function(x) { min(x$beta_hist[,2]) }, res_noisy)),
	unlist(Map(function(x) { min(x$beta_hist[,2]) }, res_noisefree))
)
y_hi = max(
	unlist(Map(function(x) { max(x$beta_hist[,2]) }, res_noisy)),
	unlist(Map(function(x) { max(x$beta_hist[,2]) }, res_noisefree))
)
y_ran = c(y_lo, y_hi)

x_ticks = round(seq(x_ran[1], x_ran[2], length.out = 10), 2)
y_ticks = round(seq(y_ran[1], y_ran[2], length.out = 10), 2)

idx = 0
for (idx_rho in seq_along(rho_levels)) {
for (idx_lambda in seq_along(lambda_levels)) {
	idx = idx + 1
	noisy_out = res_noisy[[idx]]
	noisefree_out = res_noisefree[[idx]]

	rho_y = rep(rho_levels[idx_rho], n)
	lambda_x = matrix(lambda_levels[idx_lambda], n, d1)

	printf("%d: rho_levels[%d] = %f, lambda_levels[%d] = %f\n",
		idx, idx_rho, rho_y[1], idx_lambda, lambda_x[1,1])

	g = my_contour_plot(
		noisy_out$beta_hist[,1],
		noisy_out$beta_hist[,2],
		beta_true[1], beta_true[2]) +
		xlim(x_ran[1], x_ran[2]) +
		ylim(y_ran[1], y_ran[2])
	ff = sprintf("beta_noisy_rho%03d_lambda%03d.pdf", idx_rho, idx_lambda)
	ggsave(ff, g, width = 4, height = 3)

	g = my_contour_plot(
		noisefree_out$beta_hist[,1],
		noisefree_out$beta_hist[,2],
		beta_true[1], beta_true[2]) +
		xlim(x_ran[1], x_ran[2]) +
		ylim(y_ran[1], y_ran[2])
	ff = sprintf("beta_noisefree_rho%03d_lambda%03d.pdf", idx_rho, idx_lambda)
	ggsave(ff, g, width = 4, height = 3)

	g = my_traceplot(noisy_out$sigma2_hist, "sigma2") +
		geom_hline(yintercept = sigma2_true, col = "red", linetype = 2)
	ff = sprintf("sigma2_noisy_rho%03d_lambda%03d.pdf", idx_rho, idx_lambda)
	ggsave(ff, g, width = 4, height = 3)

	g = my_traceplot(noisefree_out$sigma2_hist, "sigma2") +
		geom_hline(yintercept = sigma2_true, col = "red", linetype = 2)
	ff = sprintf("sigma2_noisefree_rho%03d_lambda%03d.pdf", idx_rho, idx_lambda)
	ggsave(ff, g, width = 4, height = 3)
}
}

# Alternative traceplots of sigma2, grouping noisy and noisefree together
R = length(res_noisy[[1]]$sigma2_hist)
dat1_plot = data.frame(x = 1:R, y = res_noisy[[1]]$sigma2_hist)
dat2_plot = data.frame(x = 1:R, y = res_noisefree[[1]]$sigma2_hist)
g = ggplot() +
	geom_line(data = dat1_plot, aes(x = x, y = y), lty = 1, col = "black") +
	geom_line(data = dat2_plot, aes(x = x, y = y), lty = 1, col = "darkgrey") +
	geom_hline(yintercept = sigma2_true, col = "red", linetype = 2) +
	xlab("Saved Draws") +
	ylab(bquote(sigma^2)) +
	theme_bw()
ff = sprintf("sigma2_rho%03d_lambda%03d.pdf", 1, 1)
ggsave(ff, g, width = 4, height = 3)

R = length(res_noisy[[2]]$sigma2_hist)
dat1_plot = data.frame(x = 1:R, y = res_noisy[[2]]$sigma2_hist)
dat2_plot = data.frame(x = 1:R, y = res_noisefree[[1]]$sigma2_hist)
g = ggplot() +
	geom_line(data = dat1_plot, aes(x = x, y = y), lty = 1, col = "black") +
	geom_line(data = dat2_plot, aes(x = x, y = y), lty = 1, col = "darkgrey") +
	geom_hline(yintercept = sigma2_true, col = "red", linetype = 2) +
	xlab("Saved Draws") +
	ylab(bquote(sigma^2)) +
	theme_bw()
ff = sprintf("sigma2_rho%03d_lambda%03d.pdf", 1, 2)
ggsave(ff, g, width = 4, height = 3)
