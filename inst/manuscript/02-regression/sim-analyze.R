library(ggplot2)
library(dplyr)
source("gibbs.R")
source("../shared/util.R")
source("../shared/plots.R")

dir.create("outputs", showWarnings = FALSE)

num_rho_levels = 3
num_sigma_levels = 3

AMSE1 = matrix(NA, num_rho_levels, num_sigma_levels)
AMSE2 = matrix(NA, num_rho_levels, num_sigma_levels)
rownames(AMSE1) = rownames(AMSE2) = sprintf("rho_%d", seq_len(num_rho_levels))
colnames(AMSE1) = colnames(AMSE2) = sprintf("sigma_%d", seq_len(num_sigma_levels))

for (idx_rho in seq_len(num_rho_levels)) {
for (idx_sigma in seq_len(num_sigma_levels)) {
	ff = sprintf("sim_rho%03d_sigma%03d/results.Rdata", idx_rho, idx_sigma)
	if (!file.exists(ff)) { next }

	out = local({
		load(ff)

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

