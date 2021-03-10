library(ggplot2)

my_traceplot = function(draws, name) {
	R = length(draws)
	dat_plot = data.frame(x = 1:R, y = draws)
	ggplot(dat_plot, aes(x = x, y = y)) +
		geom_line() +
		xlab("Saved Draws") +
		ylab(name) +
		theme_bw()
}

# Plot the posterior distribution of (beta1, beta2) to see if the relationship
# between X and Y has been recovered.
my_contour_plot = function(beta1_hat, beta2_hat, beta1_true, beta2_true) {
	dat_plot = data.frame(
		beta1 = beta1_hat,
		beta2 = beta2_hat
	)
	# pal = colorRampPalette(c("white", "red"))
	ggplot(dat_plot, aes(x = beta1, y = beta2)) +
		geom_density_2d_filled() +
		scale_fill_grey(start = 1.0, end = 0.0) +
		geom_density_2d(size = 0.25, colour = "black") +
		# scale_x_continuous(expand = c(0,0)) +
		# scale_y_continuous(expand = c(0,0)) +
		geom_vline(xintercept = beta1_true, linetype = 2) +
		geom_hline(yintercept = beta2_true, linetype = 2) +
		xlab(bquote(beta[1])) +
		ylab(bquote(beta[2])) +
		theme_bw() +
		theme(legend.position = "none",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank())
}
