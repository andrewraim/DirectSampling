# Examples of direct sampling with DGeom and Lognormal
source("functions.R")
set.seed(1234)

# ----- Example 1: Basic implementation, works okay -----
n = 1000
z = 10
mu = 3
sigma2 = 2
rho = 0.7

w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
out = direct_sampler_basic(n, w_obj, g_obj, N = 100, max_iter = 10000)

# Histogram showing number of tries needed for successful draw
dat_plot = as.data.frame(out)
g = ggplot(dat_plot, aes(x = tries)) + 
	geom_density() +
	xlab("Rejection Sampling Attempts") +
	ylab("Density") +
	theme_bw()
ggsave("example1_tries.pdf", g, width = 4, height = 3)

# Tabulate draws from the sampler
p_seq = table(out$x) / n
x_seq = as.integer(names(p_seq))

# Using draws from the sampler, evaluate the target density
# with an approximate normalizing constant
# Unnormalized target log-density function
log_f = function(x, z, mu, sigma2, rho) {
	w_obj = get_lognormal_weight(z, mu, sigma2)
	g_obj = get_dgeom_base(rho)
	w_obj$eval(x, log = TRUE) + g_obj$density(x, log = TRUE)
}
f_seq = log_f(x_seq, z, mu, sigma2, rho)

# Compare draws with density
dat_draws = data.frame(x = x_seq, p = as.numeric(p_seq))
dat_fn = data.frame(x = x_seq, p = exp(f_seq) / sum(exp(f_seq)))
g = ggplot() + 
	geom_point(data = dat_draws, aes(x = x, y = p, pch = "1")) +
	geom_line(data = dat_fn, aes(x = x, y = p)) +
	geom_point(data = dat_fn, aes(x = x, y = p, pch = "3")) +
	xlab("x") +
	ylab("Probability") +
	theme_bw() +
	scale_shape_manual(name = "Points", values = c(1, 3), breaks = c("1", "3"),
		labels = c("Empirical","Density")) +
	theme(legend.position = 'none')
ggsave("example1_draws.pdf", g, width = 4, height = 3)

# ----- Example 2: Basic implementation, works poorly -----
n = 1000
z = 244388
mu = 3.5
sigma2 = 8.5
rho = 0.7

w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
out = direct_sampler_basic(n, w_obj, g_obj, N = 100, max_iter = 10000)
table(out$tries)

# ----- Example 2: Customized implementation, works okay -----
n = 1000
z = 244388
mu = 3.5
sigma2 = 8.5
rho = 0.7

w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
x = direct_sampler(n, w_obj, g_obj, tol = 1e-8, N = 100)

# Tabulate draws from the sampler
p_seq = table(x) / n
x_seq = as.integer(names(p_seq))

# Using draws from the sampler, evaluate the target density
# with an approximate normalizing constant
log_f = function(x, z, mu, sigma2, rho) {
	w_obj = get_lognormal_weight(z, mu, sigma2)
	g_obj = get_dgeom_base(rho)
	w_obj$eval(x, log = TRUE) + g_obj$density(x, log = TRUE)
}
f_seq = log_f(x_seq, z, mu, sigma2, rho)

# Compare draws with density
dat_draws = data.frame(x = x_seq, p = as.numeric(p_seq))
dat_fn = data.frame(x = x_seq, p = exp(f_seq) / sum(exp(f_seq)))
g = ggplot() + 
	geom_point(data = dat_draws, aes(x = x, y = p, pch = "1")) +
	geom_line(data = dat_fn, aes(x = x, y = p)) +
	geom_point(data = dat_fn, aes(x = x, y = p, pch = "3")) +
	xlab("x") +
	ylab("Probability") +
	theme_bw() +
	scale_shape_manual(name = "Points", values = c(1, 3), breaks = c("1", "3"),
		labels = c("Empirical","Density")) +
	theme(legend.position = 'none')
ggsave("example2_draws.pdf", g, width = 4, height = 3)

# Find the values of u_L and u_H in this setting
step = Stepdown$new(w_obj, g_obj, tol = 1e-8, N = 100, method = "small_rects")
u = exp(step$get_log_x_vals())
printf("u_L = %g\n", u[2])
printf("u_H = %g\n", tail(u,1))

# ----- Example 3: Customized implementation, try a few more settings -----
# Get a label to use for the title of the plots
get_label = function(z, mu, sigma2, rho) {
	bquote(z == .(z) ~ "," ~ mu == .(mu) ~ "," ~
		sigma^2 == .(sigma2) ~ "," ~ rho == .(rho))
}

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; rho = 0.7
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3a_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; rho = 0.1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3b_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; rho = 0.01
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3c_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -10; mu = 3.5; sigma2 = 8.5; rho = 0.7
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3d_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -10; mu = 3.5; sigma2 = 8.5; rho = 0.1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3e_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -10; mu = 3.5; sigma2 = 8.5; rho = 0.01
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example3f_draws.pdf", gg, width = 4, height = 3)

# ----- Example 4 -----
n = 50000; z = -244; mu = 3.5; sigma2 = 8.5; rho = 0.9
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
print(gg)

n = 50000; z = -244; mu = 3.5; sigma2 = 8.5; rho = 0.1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
print(gg)

n = 50000; z = -244; mu = 3.5; sigma2 = 8.5; rho = 0.01
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
print(gg)
ggsave("example4_draws.pdf", gg, width = 4, height = 3)
