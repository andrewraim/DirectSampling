library(DirectSampling)

# Examples of direct sampling with DscNorm and Lognormal
source("functions.R")
set.seed(1234)

# Get a label to use for the title of the plots
get_label = function(z, mu, sigma2, tau) {
	bquote(z == .(z) ~ "," ~ mu == .(mu) ~ "," ~
		sigma^2 == .(sigma2) ~ "," ~ tau == .(tau))
}

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; tau = 5
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5a_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; tau = 1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5b_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = 244388; mu = 3.5; sigma2 = 8.5; tau = 0.25
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5c_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; tau = 5
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5d_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; tau = 1
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5e_draws.pdf", gg, width = 4, height = 3)

n = 50000; z = -1; mu = 3.5; sigma2 = 8.5; tau = 0.25
w_obj = get_lognormal_weight(z, mu, sigma2)
g_obj = get_dscnorm_base(tau)
gg = make_plot_discrete(n, w_obj, g_obj) +
	ggtitle(get_label(z, mu, sigma2, tau))
ggsave("example5f_draws.pdf", gg, width = 4, height = 3)
