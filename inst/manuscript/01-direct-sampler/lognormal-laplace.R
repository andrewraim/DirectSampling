# Examples of direct sampling with Laplace and Lognormal
source("functions.R")
set.seed(1234)

# Get a label to use for the title of the plots
get_label = function(z, mu, sigma2, lambda) {
	bquote(z == .(z) ~ "," ~ mu == .(mu) ~ "," ~
		sigma^2 == .(sigma2) ~ "," ~ lambda == .(lambda))
}

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; lambda = 5
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5a_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; lambda = 1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5b_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; lambda = 0.25
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5c_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; lambda = 5
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5d_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; lambda = 1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5e_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; lambda = 0.25
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example5f_draws.pdf", gg, width = 4, height = 3)
