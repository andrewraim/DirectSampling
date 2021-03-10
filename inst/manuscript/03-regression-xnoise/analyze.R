library(ggplot2)
library(dplyr)
source("gibbs.R")
source("../shared/plots.R")

dir.create("outputs", showWarnings = FALSE)

rho_levels = c(0.025, 0.10, 0.40)
sigma2_levels = c(0.25^2, 1^2, 1.25^2)

AMSE1 = matrix(NA, length(rho_levels), length(sigma2_levels))
AMSE2 = matrix(NA, length(rho_levels), length(sigma2_levels))
rownames(AMSE1) = rownames(AMSE2) = sprintf("rho_%d", seq_along(rho_levels))
colnames(AMSE1) = colnames(AMSE2) = sprintf("sigma_%d", seq_along(sigma2_levels))

for (idx_rho in seq_along(rho_levels)) {
for (idx_sigma in seq_along(sigma2_levels)) {
	ff = sprintf("sim_rho%03d_sigma%03d/results.Rdata", idx_rho, idx_sigma)
	if (!file.exists(ff)) { next }

	out = local({
		load(ff)

		# Compare quantiles of the distributions of Y
		y_hi = qlnorm(0.99, mu_true, sqrt(sigma2_true))
		plot_dat = data.frame(x = 1:n, log_y_hi = sort(log10(y_hi)))
		g = ggplot(plot_dat, aes(x = x, y = log_y_hi)) + 
			geom_line(color = "black") +
			xlab("Index (Sorted)") +
			ylab("99% Quantile of Distribution on Log10 Scale") +
			theme_bw()
		ggsave(sprintf("outputs/quantile_y_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma),
			g, height = 4, width = 4)

		mse = function(x) {
			R_keep = x$R_keep
			theta_true_mat = rep(1,R_keep) %o% c(beta_true, sigma2_true)
			theta_draws_mat = cbind(x$beta_hist, x$sigma2_hist)
			mean(rowSums((theta_draws_mat - theta_true_mat)^2))
		}

		mse1 = unlist(Map(mse, res1))
		mse2 = unlist(Map(mse, res2))

		g = ggplot() +
			geom_density(data = data.frame(x = mse1), aes(x = x), linetype = "solid") +
			geom_density(data = data.frame(x = mse2), aes(x = x), linetype = "dashed") +
			xlab("MSE") +
			ylab("Density") +
			theme_bw()
		ggsave(sprintf("outputs/mse_rho%03d_sigma%03d.pdf", idx_rho, idx_sigma),
			g, height = 4, width = 4)

		list(mse1 = mse1, mse2 = mse2)
	})
	AMSE1[idx_rho, idx_sigma] = mean(out$mse1)
	AMSE2[idx_rho, idx_sigma] = mean(out$mse2)
}
}

printf("AMSE1:\n")
print_matrix_latex(AMSE1, fmt = "%0.8f")

printf("AMSE2:\n")
print_matrix_latex(AMSE2, fmt = "%0.8f")

