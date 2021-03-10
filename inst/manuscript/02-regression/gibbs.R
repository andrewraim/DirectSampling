library(Matrix)
library(DirectSampling)

get_direct_control = function(N = 100, tol = 1e-10,
	fill_method = "small_rects")
{
	ret = list(N_direct = N, fill_method = fill_method, tol = tol)
	class(ret) = "direct_control"
	return(ret)
}

get_gibbs_control = function(save_xi = FALSE, R = 1000, burn = 0, thin = 1,
	report_period = R+1, direct = get_direct_control())
{
	ret = list(save_xi = save_xi, R = R, burn = burn, thin = thin,
		report_period = report_period, direct = direct)
	class(ret) = "gibbs_control"
	return(ret)
}

get_prior = function(a_sigma, b_sigma, a_pi, sigma2_beta, d)
{
	stopifnot(all(a_sigma > 0))
	stopifnot(all(b_sigma > 0))
	stopifnot(sigma2_beta > 0)
	ret = list(a_sigma = a_sigma, b_sigma = b_sigma, sigma2_beta = sigma2_beta)
	class(ret) = "gibbs_prior"
	return(ret)
}

get_fixed = function(beta = FALSE, sigma2 = FALSE, xi = FALSE)
{
	ret = list(beta = beta, sigma2 = sigma2, xi = xi)
	class(ret) = "gibbs_fixed"
	return(ret)
}

get_init = function(beta = NULL, sigma2 = NULL, xi = NULL, n, d)
{
	if (is.null(beta)) {
		beta = rep(0, d)
	} 
	stopifnot(length(beta) == d)

	if (is.null(sigma2)) {
		sigma2 = 1
	}
	stopifnot(length(sigma2) == 1)

	if (is.null(xi)) {
		xi = rep(0, n)
	}
	stopifnot(length(xi) == n)

	ret = list(beta = beta, sigma2 = sigma2, xi = xi)
	class(ret) = "my_init"
	return(ret)
}

my_sampler = function(y_tilde, X, rho, init, prior, control, fixed = get_fixed())
{
	n = length(y_tilde)
	d = ncol(X)

	stopifnot(length(rho) == n)

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
	xi_hist = NULL
	if (control$save_xi) {
		xi_hist = matrix(NA, R_keep, n)
	}

	# Set up initial values. xi is drawn first, so initial values for it
	# don't matter unless it is treated as fixed.
	stopifnot(class(init) == "my_init")
	beta = init$beta
	sigma2 = init$sigma2
	xi = init$xi

	# Update mu from init beta
	mu = X %*% beta

	# Set up prior
	stopifnot(class(prior) == "gibbs_prior")
	a_sigma = prior$a_sigma
	b_sigma = prior$b_sigma
	sigma2_beta = prior$sigma2_beta

	# Set up fixed parameters
	if (is.null(fixed)) { fixed = get_fixed() }
	stopifnot(class(fixed) == "gibbs_fixed")

	# Set up timers
	elapsed = list(beta = 0, sigma2 = 0, xi = 0)

	for (rep in 1:R) {
		if (rep %% report_period == 0) { logger("Starting rep %d\n", rep) }

		# Draw [xi | rest]
		st = Sys.time()
		if (!fixed$xi) {
			for (i in 1:n) {
				tryCatch({
					xi[i] = direct_sampler_lognormal_dgeom(n = 1, y_tilde[i], mu[i], sigma2,
						rho[i], tol = direct_options$tol, N = direct_options$N,
						fill_method = direct_options$fill_method)
				}, error = function(e) {
					browser()
				})
			}
		}
		elapsed$xi = elapsed$xi + as.numeric(Sys.time() - st, units = "secs")

		# Update value of [y | rest]
		if (any(xi > y_tilde)) {
			stop("Something went wrong: at least one xi > y_tilde")
		}
		log_y = log(y_tilde - xi)

		# Draw [beta | rest]
		st = Sys.time()
		if (!fixed$beta) {
			Omega = symmpart( 1/sigma2 * crossprod(X,X) + diag(1/sigma2_beta, d, d) )
			mean = solve(Omega, t(X) %*% as.matrix(log_y / sigma2))
			beta = rmvnorm_prec(1, mean, Omega)
		}
		elapsed$beta = elapsed$beta + as.numeric(Sys.time() - st, units = "secs")

		# Update mu from beta draws
		mu = X %*% beta

		# Draw [sigma2 | rest]
		st = Sys.time()
		if (!fixed$sigma2) {
			a_star = a_sigma + n/2
			b_star = b_sigma + sum((log_y - mu)^2) / 2
			sigma2 = rinvgamma(1, a_star, b_star)
		}
		elapsed$sigma2 = elapsed$sigma2 + as.numeric(Sys.time() - st, units = "secs")

		if (rep > burn && rep %% thin == 0) {
			rep_keep = rep_keep + 1
			beta_hist[rep_keep,] = beta
			sigma2_hist[rep_keep] = sigma2
			if (control$save_xi) {
				xi_hist[rep_keep,] = xi
			}
		}
	}

	ret = list(beta_hist = beta_hist, sigma2_hist = sigma2_hist,
		xi_hist = xi_hist, R_keep = R_keep, elapsed = elapsed, R = R,
		burn = burn, thin = thin)
	class(ret) = "my_fit"
	return(ret)
}

summary.my_fit = function(object, ...)
{
	pr = c(0.025, 0.5, 0.975)
	d = ncol(object$beta_hist)

	df_beta = as.data.frame(cbind(
		apply(object$beta_hist, 2, mean),
		apply(object$beta_hist, 2, sd),
		t(apply(object$beta_hist, 2, quantile, probs = pr))
	))
	rownames(df_beta) = sprintf("beta[%d]", seq_len(d))

	df_sigma2 = as.data.frame(cbind(
		mean(object$sigma2_hist),
		sd(object$sigma2_hist),
		t(quantile(object$sigma2_hist, probs = pr))
	))
	rownames(df_sigma2) = sprintf("sigma2")

	df = rbind(df_beta, df_sigma2)
	colnames(df) = c("mean", "sd", "2.5%", "50%", "97.5%")
	return(df)
}

print.my_fit = function(x, ...)
{
	cat("Summary of fit for msmix stick-breaking model\n")
	print(summary(x))

	cat("----\n")
	printf("Total iterations R: %d   burn: %d   thin: %d   Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)

	cat("----\n")
	printf("Elapsed time (Seconds):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}

