# Examples of direct sampling with DGeom and Normal
source("functions.R")
set.seed(1234)

# Get a label to use for the title of the plots
get_label = function(z, mu, sigma2, rho) {
	bquote(z == .(z) ~ "," ~ mu == .(mu) ~ "," ~
		sigma^2 == .(sigma2) ~ "," ~ lambda == .(lambda))
}

n = 50000; z = 5; mu = 0; sigma2 = 1; lambda = 1
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example7a_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 5; mu = 0; sigma2 = 1; lambda = 5
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example7b_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 5; mu = 0; sigma2 = 1; lambda = 10
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_laplace_base(lambda)
gg = make_plot_continuous(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, lambda))
ggsave("example7c_draws.pdf", gg, width = 4, height = 3)
