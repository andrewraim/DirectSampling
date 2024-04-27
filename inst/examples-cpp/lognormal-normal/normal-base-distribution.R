get_normal_base = function(sigma)
{
	# Evaluate the function
	density = function(x, log = FALSE) {
		dnorm(x, sd = sigma, log = log)
	}

	# Compute Pr(x1 < X < x2) probability where X ~ N(0,sigma^2)
	pr_interval = function(x1, x2, log = log) {
		pr = pnorm(x2, sd = sigma) - pnorm(x1, sd = sigma)
		ifelse(log, log(pr), pr)
	}

	# Quantile function of Geom truncated to (x_min, x_max)
	q_truncated = function(p, x_min = -Inf, x_max = Inf) {
		p_min = pnorm(x_min, sd = sigma)
		p_max = pnorm(x_max, sd = sigma)
		qnorm((p_max - p_min)*p + p_min, sd = sigma)
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
