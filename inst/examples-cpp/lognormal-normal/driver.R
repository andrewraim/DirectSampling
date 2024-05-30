# Direct sampling with Normal base distribution and Lognormal weight function.
# This example appears in the technical reports:
# - <https://www.census.gov/library/working-papers/2021/adrm/RRS2021-01.html>
# - <https://www.census.gov/library/working-papers/2022/adrm/RRS2022-05.html>
library(tidyverse)

set.seed(1234)

z = -10
mu = 3.5
sigma2 = 8.5
tau = 1

# ----- Sampler in C++ -----
base_dir = system.file(file.path("examples-cpp", "lognormal-normal"), package="DirectSampling")
Rcpp::sourceCpp(file.path(base_dir, "sampler.cpp"))

out1 = r_lognormal_normal(n = 50000, z, mu, sigma2, tau, N = 50,
	max_rejections = 10000)

rej1_rate = out1$rejections / (length(out1$x) + out1$rejections)

data.frame(x = out1$x) %>%
	ggplot() +
	geom_histogram(aes(x = x), bins = 25, col = "black") +
	ggtitle(sprintf("Sample after %d Rejections (%0.2f%%)",
		out1$rejections, 100 * rej1_rate)) +
	theme_minimal()

# ----- Sampler in pure R for comparison -----
library(DirectSampling)
source(file.path(base_dir, "lognormal-weight-function.R"))
source(file.path(base_dir, "normal-base-distribution.R"))

w = get_lognormal_weight(z, mu, sigma2)
g = get_normal_base(tau)
out2 = direct_sampler_ar(n = 50000, w, g, tol = 1e-8, N = 50,
	max_rejections = 10000, verbose = TRUE)

rej2_rate = out2$rejections / (length(out2$x) + out2$rejections)

data.frame(x = out2$x) %>%
	ggplot() +
	geom_histogram(aes(x = x), bins = 25, col = "black") +
	ggtitle(sprintf("Sample after %d Rejections (%0.2f%%)",
		out2$rejections, 100 * rej2_rate)) +
	theme_minimal()

# ----- Compare results -----
gg = ggplot() +
	geom_density(data = data.frame(x = out1$x), aes(x = x)) +
	geom_density(data = data.frame(x = out2$x), aes(x = x), col = "blue", lty = 2) +
	theme_linedraw()
print(gg)
