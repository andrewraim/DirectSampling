# Direct sampling with Normal base distribution and Lognormal weight function
library(tidyverse)

Rcpp::sourceCpp("sampler.cpp")

out = r_lognormal_normal(n = 10000, z = -10, mu = 3.5, sigma2 = 8.5, tau = 1,
	N = 1, max_rejections = 10000)

rej_rate = out$rejections / (length(out$x) + out$rejections)

data.frame(x = out$x) %>%
	ggplot() +
	geom_histogram(aes(x = x), bins = 25, col = "black") +
	ggtitle(sprintf("Sample after %d Rejections (%0.2f%%)",
		out$rejections, 100 * rej_rate)) +
	theme_minimal()

