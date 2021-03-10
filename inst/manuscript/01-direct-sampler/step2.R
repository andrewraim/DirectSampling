library(dplyr)
source("functions.R")

# Produce Figs 1a and 1b.
# Here u is restricted to a tiny interval.
z = 200000; mu = -2; sigma2 = 8.5; rho = 0.01

N = 5
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
step1 = Stepdown$new(w_obj, g_obj, tol = 1e-10, N = N, method = "equal_steps")
step2 = Stepdown$new(w_obj, g_obj, tol = 1e-10, N = N, method = "small_rects")

p = function(u, step, take_log = FALSE) {
	n = length(u)
	log_p = step$get_log_p
	out = numeric(n)
	for (i in 1:n) {
		out[i] = log_p(log(u[i]))
	}
	if (take_log) { return(out) } else { return(exp(out)) }
}

# ----- Produce plots -----
dat1 = data.frame(x = exp(step1$get_log_x_vals()), y = exp(step1$get_log_h_vals()))
k = nrow(dat1)
dat1_segments = cbind(dat1[-k,], dat1[-1,]) %>%
	rename(x1 = 1, y1 = 2, x2 = 3, y2 = 4)
g = ggplot() +
	geom_point(data = dat1, aes(x=x, y=y), size = 2) +
	stat_function(fun = p, colour = "red", args = list(step = step1)) +
	geom_segment(aes(x = x1, xend = x2, y = y1, yend = y1), data = dat1_segments, linetype = 2) +
	geom_segment(aes(x = x2, xend = x2, y = y1, yend = y2), data = dat1_segments, linetype = 2) +
	xlab("u") +
	ylab(bquote(P(A[u]))) +
	theme_bw()
ggsave("step2_1.pdf", g, width = 4, height = 3)

dat2 = data.frame(x = exp(step2$get_log_x_vals()), y = exp(step2$get_log_h_vals()))
k = nrow(dat2)
dat2_segments = cbind(dat2[-k,], dat2[-1,]) %>%
	rename(x1 = 1, y1 = 2, x2 = 3, y2 = 4)
g = ggplot() +
	geom_point(data = dat2, aes(x=x, y=y), size = 2) +
	stat_function(fun = p, colour = "red", args = list(step = step2)) +
	geom_segment(aes(x = x1, xend = x2, y = y1, yend = y1), data = dat2_segments, linetype = 2) +
	geom_segment(aes(x = x2, xend = x2, y = y1, yend = y2), data = dat2_segments, linetype = 2) +
	xlab("u") +
	ylab(bquote(P(A[u]))) +
	theme_bw()
ggsave("step2_2.pdf", g, width = 4, height = 3)
