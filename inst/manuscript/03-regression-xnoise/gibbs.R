library(Matrix)
library(DirectSampling)

get_direct_control = function(N = 100, tol = 1e-10, fill_method = "small_rects")
{
	ret = list(N_direct = N, fill_method = fill_method, tol = tol)
	class(ret) = "direct_control"
	return(ret)
}

get_gibbs_control = function(save_xi_y = FALSE, save_xi_x = FALSE,
	R = 1000, tol_beta = 1e-6, burn = 0, thin = 1, report_period = R+1,
	direct_y = get_direct_control(), direct_x = get_direct_control())
{
	ret = list(save_xi_y = save_xi_y, save_xi_x = save_xi_x, R = R,
		burn = burn, thin = thin, tol_beta = tol_beta,
		report_period = report_period, direct_y = direct_y, direct_x = direct_x)
	class(ret) = "gibbs_control"
	return(ret)
}

get_prior = function(a_sigma, b_sigma, sigma2_beta, d)
{
	stopifnot(all(a_sigma > 0))
	stopifnot(all(b_sigma > 0))
	stopifnot(sigma2_beta > 0)
	ret = list(a_sigma = a_sigma, b_sigma = b_sigma, sigma2_beta = sigma2_beta)
	class(ret) = "gibbs_prior"
	return(ret)
}

get_fixed = function(beta = FALSE, sigma2 = FALSE, xi_y = FALSE, xi_x = FALSE)
{
	ret = list(beta = beta, sigma2 = sigma2, xi_y = xi_y, xi_x = xi_x)
	class(ret) = "gibbs_fixed"
	return(ret)
}

get_init = function(beta = NULL, sigma2 = NULL, xi_y = NULL, xi_x = NULL, n, d)
{
	if (is.null(beta)) {
		beta = rep(0, d)
	} 
	stopifnot(length(beta) == d)

	if (is.null(sigma2)) {
		sigma2 = 1
	}
	stopifnot(length(sigma2) == 1)

	if (is.null(xi_y)) {
		xi_y = rep(0, n)
	}
	stopifnot(length(xi_y) == n)

	if (is.null(xi_x)) {
		xi_x = matrix(0, n, d)
	}
	stopifnot(dim(xi_x) == c(n,d))

	ret = list(beta = beta, sigma2 = sigma2, xi_y = xi_y, xi_x = xi_x)
	class(ret) = "gibbs_init"
	return(ret)
}

my_sampler = function(y_tilde, X_tilde, rho_y, lambda_x, init, prior, control,
	fixed = get_fixed())
{
	n = length(y_tilde)
	d = ncol(X_tilde)
	d1 = ncol(lambda_x)

	stopifnot(length(rho_y) == n)
	stopifnot(dim(lambda_x) == c(n,d1))

	stopifnot(class(control) == "gibbs_control")
	R = control$R
	burn = control$burn
	thin = control$thin
	report_period = control$report_period
	tol_beta = control$tol_beta
	direct_options_y = control$direct_y
	direct_options_x = control$direct_x
	
	stopifnot(class(direct_options_y) == "direct_control")
	stopifnot(class(direct_options_x) == "direct_control")
	stopifnot(tol_beta > 0)

	rep_keep = 0
	R_keep = ceiling((R - burn) / thin)

	# Set up histories
	beta_hist = matrix(NA, R_keep, d)
	sigma2_hist = numeric(R_keep)
	xi_y_hist = NULL
	xi_x_hist = NULL
	if (control$save_xi_y) {
		xi_y_hist = matrix(NA, R_keep, n)
	}
	if (control$save_xi_x) {
		xi_x_hist = array(NA, dim = c(R_keep, n, d1))
	}

	# Set up initial values. xi_y and xi_x are drawn first, so initial values
	# for them don't matter unless they are treated as fixed
	stopifnot(class(init) == "gibbs_init")
	beta = init$beta
	sigma2 = init$sigma2
	xi_y = init$xi_y
	xi_x = init$xi_x

	# Update value of X
	X = X_tilde - xi_x

	# Update mu from initial beta
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
	elapsed = list(beta = 0, sigma2 = 0, xi_y = 0, xi_x = 0)

	for (rep in 1:R) {
		if (rep %% report_period == 0) { logger("Starting rep %d\n", rep) }

		# Draw [xi_y | rest]
		st = Sys.time()
		if (!fixed$xi_y) {
			for (i in 1:n) {
				xi_y[i] = direct_sampler_lognormal_dgeom(n = 1, y_tilde[i], mu[i], sigma2,
					rho_y[i], N = direct_options_y$N, tol = direct_options_y$tol,
					fill_method = direct_options_y$fill_method)
			}
		}
		elapsed$xi_y = elapsed$xi_y + as.numeric(Sys.time() - st, units = "secs")

		# Update value of [y | rest]
		if (any(xi_y > y_tilde)) {
			stop("Something went wrong: at least one xi_y > y_tilde")
		}
		log_y = log(y_tilde - xi_y)

		# Draw [xi_x | rest]
		st = Sys.time()
		if (!fixed$xi_x) {
			for (j in 1:d1) {
				if (abs(beta[j]) > tol_beta) {
					idx_minus_j = setdiff(1:d, j)
					z_j = log_y - X[,idx_minus_j,drop = FALSE] %*% beta[idx_minus_j]
					mean_j = z_j / beta[j]
					tau2_j = sigma2 / beta[j]^2
					for (i in 1:n) {
						xi_x[i,j] = direct_sampler_normal_laplace(n = 1, X_tilde[i,j],
							mean_j[i], tau2_j, lambda_x[i,j], N = direct_options_x$N,
							tol = direct_options_x$tol,
							fill_method = direct_options_x$fill_method)
					}
				} else {
					for (i in 1:n) {
						xi_x[i,j] = r_dgeom(n = 1, lambda_x[i,j])
					}
				}

				# Update value of X[,j]
				X[,j] = X_tilde[,j] - xi_x[,j]
			}
		}
		elapsed$xi_x = elapsed$xi_x + as.numeric(Sys.time() - st, units = "secs")

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
			if (control$save_xi_y) {
				xi_y_hist[rep_keep,] = xi_y
			}
			if (control$save_xi_x) {
				xi_x_hist[rep_keep,,] = xi_x[,1:d1]
			}
		}
	}

	ret = list(beta_hist = beta_hist, sigma2_hist = sigma2_hist,
		xi_y_hist = xi_y_hist, xi_x_hist = xi_x_hist, R_keep = R_keep,
		elapsed = elapsed, R = R, burn = burn, thin = thin)
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
	cat("Summary of fit\n")
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
