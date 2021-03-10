# Examples of direct sampling with DGeom and Normal
source("functions.R")
set.seed(1234)

# Get a label to use for the title of the plots
get_label = function(z, mu, sigma2, rho) {
	bquote(z == .(z) ~ "," ~ mu == .(mu) ~ "," ~
		sigma^2 == .(sigma2) ~ "," ~ rho == .(rho))
}

n = 50000; z = 5; mu = 0; sigma2 = 1; rho = 0.9
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example6a_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 5; mu = 0; sigma2 = 1; rho = 0.2
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example6b_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 5; mu = 0; sigma2 = 1; rho = 0.01
w_obj = get_normal_weight(z, mu, sigma2)
g_obj = get_dgeom_base(rho)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, rho))
ggsave("example6c_draws.pdf", gg, width = 4, height = 3)

