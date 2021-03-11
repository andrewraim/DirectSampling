library(Rcpp)
library(ggplot2)
library(gridExtra)

# ----- Setup for customized direct sampler -----
# The customized version of the direct sampler is written more generally. It
# takes a "weight" object and a "base" object as input. The following functions
# produce these kinds of objects.

get_lognormal_weight = function(z, mu, sigma2)
{
	# log_c is the maximum value of the function log[w(x)]
	log_c = -(mu - sigma2 / 2)

	# Evaluate the function
	eval = function(x, log = FALSE) {
		n = length(x)
		out = rep(-Inf, n)
		idx = which(x < z)
		out[idx] = -log(z-x[idx]) - (log(z-x[idx]) - mu)^2 / (2*sigma2)
		if (log) { return(out) } else { return(exp(out)) }
	}

	# Return the roots of the equation w(x) = a, which is equivalent to
	# log(w(x)) = log(a). Roots are returned in increasing order.
	roots = function(log_a) {
		if (sigma2 < 2*(mu + log_a) && log_a <= log_c) {
			# Handle this case to prevent NaN issues due to numerical error
			x1 = z - exp(mu - sigma2)
			x2 = z - exp(mu - sigma2)
		} else {
			x1 = z - exp((mu - sigma2) + sqrt(sigma2 * (sigma2 - 2*(mu + log_a))))
			x2 = z - exp((mu - sigma2) - sqrt(sigma2 * (sigma2 - 2*(mu + log_a))))
		}

		# Edge case: If z is an integer and the larger root is numerically close to
		# z, make the root smaller by an epsilon. This helps to ensure that the
		# endpoint z is not part of the support of p(u).
		if (z - floor(z) < 1e-10 && z - x2 < 1e-6) {
			x2 = z - 1e-6
		}

		c(x1, x2)
	}

	ret = list(log_c = log_c, roots = roots, eval = eval)
	class(ret) = "weight"
	return(ret)
}

get_normal_weight = function(z, mu, sigma2)
{
	# log_c is the maximum value of the function log[w(x)]
	log_c = 0

	# Evaluate the function
	eval = function(x, log = FALSE) {
		out = -1/(2*sigma2) * (z - x - mu)^2
		if (log) { return(out) } else { return(exp(out)) }
	}

	# Return the roots of the equation w(x) = a, which is equivalent to
	# log(w(x)) = log(a). Roots are returned in increasing order.
	roots = function(log_a) {
		x1 = z - (mu + sqrt(-2 * sigma2 * log_a))
		x2 = z - (mu - sqrt(-2 * sigma2 * log_a))
		c(x1, x2)
	}

	ret = list(log_c = log_c, roots = roots, eval = eval)
	class(ret) = "weight"
	return(ret)
}

get_normal_base = function(sigma)
{
	# Evaluate the function
	density = function(x, log = FALSE) {
		dnorm(x, 0, sigma, log)
	}

	# Compute Pr(x1 < X <= x2) probability where X ~ N(0, sigma2)
	pr_interval = function(x1, x2) {
		pnorm(x2, 0, sigma) - pnorm(x1, 0, sigma)
	}

	# Quantile function of Normal truncated to [x_min, x_max]
	# As with dgeom, infinite x is handled before returning.
	q_truncated = function(p, x_min = -Inf, x_max = Inf) {
		p_min = pnorm(x_min, 0, sigma)
		p_max = pnorm(x_max, 0, sigma)
		x = qnorm((p_max - p_min)*p + p_min, 0, sigma)
		max(x_min, min(x, x_max))
	}

	r_truncated = function(n, x_min = -Inf, x_max = Inf) {
		u = runif(n)
		x = numeric(n)
		for (i in 1:n) {
			x[i] = q_truncated(u[i], x_min, x_max)
		}
		return(x)
	}

	ret = list(pr_interval = pr_interval, q_truncated = q_truncated,
		r_truncated = r_truncated, density = density)
	class(ret) = "base"
	return(ret)
}

get_dgeom_base = function(rho)
{
	# Evaluate the function
	density = function(x, log = FALSE) {
		d_dgeom(x, rho, log)
	}

	# Compute Pr(x1 < X < x2) probability where X ~ DGeom(rho)
	# This calculation should exclude endpoints if they are integers
	pr_interval = function(x1, x2) {
		a = floor(x1 + 1)
		b = ceiling(x2 - 1)
		p_dgeom(b, rho) - p_dgeom(a - 1, rho)
	}

	# Quantile function of DGeom truncated to (x_min, x_max)
	#
	# The intermediate result x can be infinite; for example, suppose
	# x_min and x_max are both large enough so that Pr(x >= x_min) and
	# Pr(x >= x_max) are numerically equal to 1. Then we return
	# q_dgeom(1, rho) = Inf.
	q_truncated = function(p, x_min = -Inf, x_max = Inf) {
		a = floor(x_min + 1)
		b = ceiling(x_max - 1)
		p_min = p_dgeom(a - 1, rho)
		p_max = p_dgeom(b, rho)
		x = q_dgeom((p_max - p_min)*p + p_min, rho)
		max(ceiling(x_min), min(x, floor(x_max)))
	}

	r_truncated = function(n, x_min = -Inf, x_max = Inf) {
		u = runif(n)
		x = numeric(n)
		for (i in 1:n) {
			x[i] = q_truncated(u[i], x_min, x_max)
		}
		return(x)
	}
	
	ret = list(pr_interval = pr_interval, q_truncated = q_truncated,
		r_truncated = r_truncated, density = density)
	class(ret) = "base"
	return(ret)
}

get_dscnorm_base = function(tau, tol = 1e-10)
{
	# Evaluate the function
	density = function(x, log = FALSE) {
		d_dscnorm(x, tau, tol = tol, take_log = log)
	}

	# Compute Pr(x1 < X <= x2) probability where X ~ DscNorm(tau)
	# This calculation should exclude endpoints if they are integers
	pr_interval = function(x1, x2) {
		a = floor(x1 + 1)
		b = ceiling(x2 - 1)
		p_dscnorm(b, tau, tol) - p_dscnorm(a - 1, tau, tol)
	}

	# Quantile function of DscNorm truncated to (x_min, x_max)
	q_truncated = function(p, x_min = -Inf, x_max = Inf) {
		a = floor(x_min + 1)
		b = ceiling(x_max - 1)
		p_min = p_dscnorm(a - 1, tau, tol)
		p_max = p_dscnorm(b, tau, tol)
		x = q_dscnorm((p_max - p_min)*p + p_min, tau, tol)
		max(ceiling(x_min), min(x, floor(x_max)))
	}

	r_truncated = function(n, x_min = -Inf, x_max = Inf) {
		u = runif(n)
		x = numeric(n)
		for (i in 1:n) {
			x[i] = q_truncated(u[i], x_min, x_max)
		}
		return(x)
	}
	
	ret = list(pr_interval = pr_interval, q_truncated = q_truncated,
		r_truncated = r_truncated, density = density)
	class(ret) = "base"
	return(ret)
}

get_laplace_base = function(lambda)
{
	# Evaluate the function
	density = function(x, log = FALSE) {
		d_laplace(x, 0, lambda, log)
	}

	# Compute Pr(x1 < X < x2) probability where X ~ Laplace(0, lambda)
	pr_interval = function(x1, x2) {
		p_laplace(x2, 0, lambda) - p_laplace(x1, 0, lambda)
	}

	# Quantile function of Laplace truncated to (x_min, x_max)
	# As with dgeom, infinite x is handled before returning.
	q_truncated = function(p, x_min = -Inf, x_max = Inf) {
		p_min = p_laplace(x_min, 0, lambda)
		p_max = p_laplace(x_max, 0, lambda)
		x = q_laplace((p_max - p_min)*p + p_min, 0, lambda)
		max(x_min, min(x, x_max))
	}

	r_truncated = function(n, x_min = -Inf, x_max = Inf) {
		u = runif(n)
		x = numeric(n)
		for (i in 1:n) {
			x[i] = q_truncated(u[i], x_min, x_max)
		}
		return(x)
	}

	ret = list(pr_interval = pr_interval, q_truncated = q_truncated,
		r_truncated = r_truncated, density = density)
	class(ret) = "base"
	return(ret)
}

# ----- Function to draw a sample and make a plot -----
make_plot_discrete = function(n, w_obj, g_obj, tol = 1e-8, N = 100)
{
	x = direct_sampler(n, w_obj, g_obj, tol, N)

	# Tabulate draws from the sampler
	p_seq = table(x) / n
	x_seq = as.integer(names(p_seq))

	# Using draws from the sampler, evaluate the target density
	# with an approximate normalizing constant
	f_seq = w_obj$eval(x_seq, log = TRUE) + g_obj$density(x_seq, log = TRUE)

	# Compare draws with density
	dat_draws = data.frame(x = x_seq, p = as.numeric(p_seq))
	dat_fn = data.frame(x = x_seq, p = exp(f_seq) / sum(exp(f_seq)))
	ggplot() +
		geom_point(data = dat_draws, aes(x = x, y = p, pch = "1")) +
		geom_line(data = dat_fn, aes(x = x, y = p)) +
		geom_point(data = dat_fn, aes(x = x, y = p, pch = "3")) +
		xlab("x") +
		ylab("Probability") +
		theme_bw() +
		scale_shape_manual(name = "Points", values = c(1, 3), breaks = c("1", "3"),
			labels = c("Empirical","Density")) +
		theme(legend.position = 'none',
			  plot.title = element_text(size = rel(1)))
}

make_plot_continuous = function(n, w_obj, g_obj, tol = 1e-8, N = 100)
{
	x = direct_sampler(n, w_obj, g_obj, tol, N)

	# Integrate to obtain the normalizing constant
	fn = function(x) { w_obj$eval(x) * g_obj$density(x) }
	norm_const = integrate(fn, lower = -Inf, upper = Inf)$value
	f_normalized = function(x) { w_obj$eval(x) * g_obj$density(x) / norm_const }

	# Plot the samples vs the density
	ggplot() +
		geom_density(data = data.frame(x = x), aes(x = x)) +
		stat_function(fun = f_normalized, linetype = 2) +
		xlab("x") +
		ylab("Probability") +
		theme_bw() +
		theme(legend.position = 'none',
			  plot.title = element_text(size = rel(1)))
}
