source("../shared/util.R")
source("../shared/invgamma.R")
source("../shared/mvnorm.R")
Rcpp::sourceCpp("sampler.cpp")

get_direct_control = function(N = 100, tol = 1e-10, fill_method = "small_rects",
	max_rejections = 0, priority_weight = 0.5)
{
	ret = list(N_direct = N, fill_method = fill_method, tol = tol,
		max_rejections = max_rejections, priority_weight = priority_weight)
	class(ret) = "direct_control"
	return(ret)
}

get_gibbs_control = function(save_s = FALSE, R = 1000, burn = 0, thin = 1,
	report_period = R+1, direct = get_direct_control())
{
	ret = list(save_s = save_s, R = R, burn = burn, thin = thin,
		report_period = report_period, direct = direct)
	class(ret) = "gibbs_control"
	return(ret)
}

get_prior = function(nu_min, nu_max, a_sigma, b_sigma, sig2beta)
{
	stopifnot(nu_min > 0)
	stopifnot(nu_max > nu_min)
	ret = list(nu_min = nu_min, nu_max = nu_max, a_sigma = a_sigma,
		b_sigma = b_sigma, sig2beta = sig2beta)
	class(ret) = "gibbs_prior"
	return(ret)
}

get_fixed = function(beta = FALSE, sigma2 = FALSE, nu = FALSE, s = FALSE)
{
	ret = list(beta = beta, sigma2 = sigma2, nu = nu, s = s)
	class(ret) = "gibbs_fixed"
	return(ret)
}

get_init = function(n = n, d = d, beta = NULL, sigma2 = NULL, nu = NULL, s = NULL)
{
	if (is.null(beta)) {
		beta = numeric(d)
	} 
	stopifnot(length(beta) == d)

	if (is.null(sigma2)) {
		sigma2 = 1
	}
	stopifnot(length(sigma2) == 1)

	if (is.null(nu)) {
		nu = 1
	}
	stopifnot(length(nu) == 1)

	if (is.null(s)) {
		s = numeric(n)
	}
	stopifnot(length(s) == n)

	ret = list(beta = beta, sigma2 = sigma2, nu = nu, s = s)
	class(ret) = "gibbs_init"
	return(ret)
}

gibbs_sampler = function(y, X, init, prior, control, fixed = get_fixed())
{
	n = length(y)
	d = ncol(X)
	stopifnot(n == nrow(X))
	
	stopifnot(class(control) == "gibbs_control")
	R = control$R
	burn = control$burn
	thin = control$thin
	report_period = control$report_period
	direct_options = control$direct
	stopifnot(class(direct_options) == "direct_control")

	rep_keep = 0
	R_keep = ceiling((R - burn) / thin)

	# Set up histories
	beta_hist = matrix(NA, R_keep, d)
	sigma2_hist = numeric(R_keep)
	nu_hist = numeric(R_keep)
	nu_rejections_hist = numeric(R_keep)

	s_hist = NULL
	if (control$save_s) {
		s_hist = matrix(NA, R_keep, n)
	}

	# Set up initial values. xi1 is drawn first, so initial values for them
	# don't matter unless they are treated as fixed.
	stopifnot(class(init) == "gibbs_init")
	beta = init$beta
	sigma2 = init$sigma2
	nu = init$nu
	nu_rejections = 0
	s = rep(NA, n)

	mu = as.numeric(X %*% beta)

	# Set up prior
	stopifnot(class(prior) == "gibbs_prior")
	stopifnot(!is.null(prior$nu_min))
	stopifnot(!is.null(prior$nu_max))
	stopifnot(!is.null(prior$a_sigma))
	stopifnot(!is.null(prior$b_sigma))
	stopifnot(!is.null(prior$sig2beta))

	# Set up fixed parameters
	if (is.null(fixed)) { fixed = get_fixed() }
	stopifnot(class(fixed) == "gibbs_fixed")

	# Set up timers
	elapsed = list(beta = 0, sigma2 = 0, nu = 0, s = 0)

	for (rep in 1:R) {
		if (rep %% report_period == 0) { logger("Starting rep %d\n", rep) }

		# Draw [ s | rest ]
		st = Sys.time()
		if (!fixed$s) {
			a_star = (nu + 1) / 2
			b_star = (y - mu)^2 / 2 + nu * sigma2 / 2
			s = 1 / rgamma(n, a_star, b_star)
		}
		elapsed$s = elapsed$s + as.numeric(Sys.time() - st, units = "secs")

		# Draw [ beta | rest ]
		st = Sys.time()
		if (!fixed$beta) {
			Omega = crossprod(X, 1/s * X) + diag(x = 1/prior$sig2beta, nrow = d)
			mean = as.numeric( solve(Omega, crossprod(X, y/s)) )
			beta = rmvnorm_prec(1, mean, Omega)
			mu = as.numeric(X %*% beta)
		}
		elapsed$beta = elapsed$beta + as.numeric(Sys.time() - st, units = "secs")

		# Draw [ sigma2 | rest ]
		st = Sys.time()
		if (!fixed$sigma2) {
			a_star = prior$a_sigma + n * nu / 2
			b_star = prior$b_sigma + nu / 2 * sum(1/s)
			sigma2 = rgamma(1, a_star, b_star)
		}
		elapsed$sigma2 = elapsed$sigma2 + as.numeric(Sys.time() - st, units = "secs")

		# Draw [ nu | rest ]
		st = Sys.time()
		if (!fixed$nu) {
			A = 1/2 * sum(log(s / sigma2)) + 1/2 * sum(sigma2 / s)

			# C++ version of the sampler
			ds_out = r_tdist_unif(n = 1, m = n, A = A,
				nu_min = prior$nu_min, nu_max = prior$nu_max,
				tol = direct_options$tol, N = direct_options$N,
				max_rejections = direct_options$max_rejections)
			nu = ds_out$x
			nu_rejections = ds_out$rejections
		}
		elapsed$nu = elapsed$nu + as.numeric(Sys.time() - st, units = "secs")

		if (rep > burn && rep %% thin == 0) {
			rep_keep = rep_keep + 1
			beta_hist[rep_keep,] = beta
			sigma2_hist[rep_keep] = sigma2
			nu_hist[rep_keep] = nu
			nu_rejections_hist[rep_keep] = nu_rejections
			if (control$save_s) {
				s_hist[rep_keep,] = s
			}
		}
	}

	ret = list(beta_hist = beta_hist, sigma2_hist = sigma2_hist, nu_hist = nu_hist,
		s_hist = s_hist, R_keep = R_keep, elapsed = elapsed, R = R, burn = burn,
		thin = thin, nu_rejections_hist = nu_rejections_hist)
	class(ret) = "gibbs_fit"
	return(ret)
}

summary.gibbs_fit = function(object, ...)
{
	pr = c(0.025, 0.5, 0.975)

	df_beta = as.data.frame(cbind(
		apply(object$beta_hist, 2, mean),
		apply(object$beta_hist, 2, sd),
		t(apply(object$beta_hist, 2, quantile, probs = pr))
	))
	d = nrow(df_beta)

	df_sigma2 = as.data.frame(cbind(
		mean(object$sigma2_hist),
		sd(object$sigma2_hist),
		t(quantile(object$sigma2_hist, probs = pr))
	))

	df_nu = as.data.frame(cbind(
		mean(object$nu_hist),
		sd(object$nu_hist),
		t(quantile(object$nu_hist, probs = pr))
	))

	df = rbind(df_beta, df_sigma2, df_nu)
	rownames(df) = c(sprintf("beta%d", 1:d), "sigma2", "nu")
	colnames(df) = c("mean", "sd", "2.5%", "50%", "97.5%")
	return(df)
}

print.gibbs_fit = function(x, ...)
{
	cat("Summary of fit\n")
	print(summary(x))

	cat("----\n")
	printf("Total iterations R: %d   burn: %d   thin: %d   Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)
	
	printf("Rejections in nu step: %d  Total proposals: %d  Rejection rate: %0.4f%%\n",
		sum(x$nu_rejections_hist),
		sum(x$nu_rejections_hist) + x$R,
		100 * sum(x$nu_rejections_hist) / (sum(x$nu_rejections_hist) + x$R))

	cat("----\n")
	printf("Elapsed time (Seconds):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}
