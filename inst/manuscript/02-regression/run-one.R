library(ggplot2)
source("gibbs.R")
source("../shared/plots.R")

set.seed(1234)

n = 200
sigma2_levels = c(0.25^2, 5^2)
rho_levels = c(0.01, 0.1)
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

# ----- Generate some covariates -----
d = length(beta_true)
X = matrix(NA, n, d)
for (l in 1:d) {
	X[,l] = rnorm(n)
}

res_noisy = list()
res_noisefree = list()
idx = 0

# ----- Fit one simulated dataset at each setting -----
for (idx_rho in seq_along(rho_levels)) {
for (idx_sigma in seq_along(sigma2_levels)) {
	sigma2_true = sigma2_levels[idx_sigma]
	rho = rep(rho_levels[idx_rho], n)

	#  Generate data from the assumed model
	mu_true = X %*% beta_true
	y_true = rlnorm(n, mu_true, sqrt(sigma2_true))
	xi_true = rgeom(n, rho) - rgeom(n, rho)
	y_tilde = y_true + xi_true

	# Setup for both samplers
	prior = get_prior(a_sigma = 2, b_sigma = 10, sigma2_beta = 10^2, d)
	control = get_gibbs_control(save_xi = TRUE, R = 2000, burn = 1000,
		thin = 1, report_period = 500, direct = get_direct_control(N = 100))

	logger("Fitting noisy obs\n")
	init = get_init(n = n, d = d)
	sampler_out = my_sampler(y_tilde, X, rho, init, prior, control)

	logger("Fitting noiseless obs\n")
	init = get_init(xi = rep(0,n), n = n, d = d)
	fixed = get_fixed(xi = TRUE)
	noisefree_out = my_sampler(y_true, X, rho = rep(1,n), init, prior,
		control, fixed)

	idx = idx + 1
	res_noisy[[idx]] = sampler_out
	res_noisefree[[idx]] = noisefree_out
}
}

# ----- Plot results -----
idx = 0
for (idx_rho in seq_along(rho_levels)) {
for (idx_sigma in seq_along(sigma2_levels)) {
	idx = idx + 1
	noisy_out = res_noisy[[idx]]
	noisefree_out = res_noisefree[[idx]]

	sigma2_true = sigma2_levels[idx_sigma]
	rho = rep(rho_levels[idx_rho], n)
	
	printf("%d: rho_levels[%d] = %f, sigma2_levels[%d] = %f\n",
		idx, idx_rho, rho[1], idx_sigma, sigma2_true)

	# Common x- and y-limits for the contour plots
	x_ran = range(noisy_out$beta_hist[,1], noisefree_out$beta_hist[,1])
	y_ran = range(noisy_out$beta_hist[,2], noisefree_out$beta_hist[,2])
	x_ticks = round(seq(x_ran[1], x_ran[2], length.out = 10), 2)
	y_ticks = round(seq(y_ran[1], y_ran[2], length.out = 10), 2)

	g = my_contour_plot(
		noisy_out$beta_hist[,1],
		noisy_out$beta_hist[,2],
		beta_true[1], beta_true[2]) +
		xlim(x_ran[1], x_ran[2]) +
		ylim(y_ran[1], y_ran[2])
	ff = sprintf("beta_noisy_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma)
	ggsave(ff, g, width = 4, height = 3)

	g = my_contour_plot(
		noisefree_out$beta_hist[,1],
		noisefree_out$beta_hist[,2],
		beta_true[1], beta_true[2]) +
		xlim(x_ran[1], x_ran[2]) +
		ylim(y_ran[1], y_ran[2])
	ff = sprintf("beta_noisefree_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma)
	ggsave(ff, g, width = 4, height = 3)

	g = my_traceplot(noisy_out$sigma2_hist, "sigma2") +
		geom_hline(yintercept = sigma2_true, col = "red", linetype = 2)
	ff = sprintf("sigma2_noisy_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma)
	ggsave(ff, g, width = 4, height = 3)

	g = my_traceplot(noisefree_out$sigma2_hist, "sigma2") +
		geom_hline(yintercept = sigma2_true, col = "red", linetype = 2)
	ff = sprintf("sigma2_noisefree_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma)
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
	geom_hline(yintercept = sigma2_levels[1], col = "red", linetype = 2) +
	xlab("Saved Draws") +
	ylab(bquote(sigma^2)) +
	theme_bw()
ff = sprintf("sigma2_rho%03d_sigma%03d.pdf", 1, 1)
ggsave(ff, g, width = 4, height = 3)

R = length(res_noisy[[2]]$sigma2_hist)
dat1_plot = data.frame(x = 1:R, y = res_noisy[[2]]$sigma2_hist)
dat2_plot = data.frame(x = 1:R, y = res_noisefree[[2]]$sigma2_hist)
g = ggplot() +
	geom_line(data = dat1_plot, aes(x = x, y = y), lty = 1, col = "black") +
	geom_line(data = dat2_plot, aes(x = x, y = y), lty = 1, col = "darkgrey") +
	geom_hline(yintercept = sigma2_levels[2], col = "red", linetype = 2) +
	xlab("Saved Draws") +
	ylab(bquote(sigma^2)) +
	theme_bw()
ff = sprintf("sigma2_rho%03d_sigma%03d.pdf", 1, 2)
ggsave(ff, g, width = 4, height = 3)
