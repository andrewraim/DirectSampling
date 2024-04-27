library(DirectSampling)
library(Matrix)

source("../shared/util.R")
source("../shared/invgamma.R")
source("../shared/mvnorm.R")
Rcpp::sourceCpp("sampler.cpp")

get_direct_control = function(N = 100, max_rejections = 0)
{
	ret = list(N = N, max_rejections = max_rejections)
	class(ret) = "direct_control"
	return(ret)
}

get_gibbs_control = function(save_eta = FALSE, R = 1000, burn = 0, thin = 1,
	report_period = R+1, direct = get_direct_control())
{
	ret = list(save_eta = save_eta, R = R, burn = burn, thin = thin,
		report_period = report_period, direct = direct)
	class(ret) = "gibbs_control"
	return(ret)
}

get_prior = function(M_sigma, M_tau, a_rho, b_rho, sigma2_beta, d)
{
	stopifnot(all(M_sigma > 0))
	stopifnot(all(M_tau > 0))
	stopifnot(all(a_rho > 0))
	stopifnot(all(b_rho > 0))
	stopifnot(sigma2_beta > 0)
	ret = list(M_sigma = M_sigma, M_tau = M_tau, a_rho = a_rho, b_rho = b_rho,
		sigma2_beta = sigma2_beta)
	class(ret) = "gibbs_prior"
	return(ret)
}

get_fixed = function(beta = FALSE, eta = FALSE, sigma2 = FALSE, tau2 = FALSE, rho = FALSE)
{
	ret = list(beta = beta, eta = eta, sigma2 = sigma2, tau2 = tau2, rho = rho)
	class(ret) = "gibbs_fixed"
	return(ret)
}

get_init = function(beta = NULL, eta = NULL, sigma2 = NULL, tau2 = NULL, rho = NULL, d, k)
{
	if (is.null(beta)) {
		beta = rep(0, d)
	} 
	stopifnot(length(beta) == d)

	if (is.null(eta)) {
		eta = rep(0, k)
	} 
	stopifnot(length(eta) == k)

	if (is.null(sigma2)) {
		sigma2 = 1
	}
	stopifnot(length(sigma2) == 1)

	if (is.null(tau2)) {
		tau2 = 1
	}
	stopifnot(length(tau2) == 1)

	if (is.null(rho)) {
		rho = 0.9
	}
	stopifnot(length(rho) == 1)

	ret = list(beta = beta, eta = eta, sigma2 = sigma2, tau2 = tau2, rho = rho)
	class(ret) = "gibbs_init"
	return(ret)
}

gibbs_sampler = function(y, X, S, A, init, prior, control, fixed = get_fixed())
{
	n = length(y)
	d = ncol(X)
	k = ncol(S)

	stopifnot(nrow(X) == n)
	stopifnot(nrow(S) == n)
	stopifnot(nrow(A) == k)
	stopifnot(ncol(A) == k)

	# Precompute these matrices to avoid repeated computation
	D = Diagonal(x = rowSums(A))

	lambda = local({
		eig = eigen(A)
		U = eig$vectors
		phi = eig$values
		phi_half = sqrt(as.complex(phi))
		V = as.matrix( t(U) %*% solve(D, U) )
		Q = Re(phi_half * (V * phi_half))
		lambda = Re(eigen(Q)$values)
		sort(lambda, decreasing = TRUE)
	})

	# The following should be equal:
	# - det(D - rho*A)
	# - prod(diag(D)) * det(diag(k) - rho*Q)
	# - prod(diag(D)) * prod(1 - rho * lambda)

	# An alternate lambda computation based on QR decomposition. This gives the
	# same overall determinant, but the eigenvalues can be complex. Therefore,
	# the one above appears simpler to use.
	# qr_out = qr(as.matrix(A))
	# Q = qr.R(qr_out) %*% solve(D, qr.Q(qr_out))
	# lambda = eigen(Q)$values

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
	tau2_hist = numeric(R_keep)
	rho_hist = numeric(R_keep)
	rho_rejections_hist = numeric(R_keep)
	eta_hist = NULL
	if (control$save_eta) {
		eta_hist = matrix(NA, R_keep, k)
	}

	# Set up initial values
	stopifnot(class(init) == "gibbs_init")
	beta = init$beta
	eta = init$eta
	sigma2 = init$sigma2
	tau2 = init$tau2
	rho = init$rho
	rho_rejections = 0

	# Update quantities based on init values
	mu = X %*% beta
	phi = S %*% eta
	H = as.matrix(D - rho*A)

	# Set up prior
	stopifnot(class(prior) == "gibbs_prior")
	M_sigma = prior$M_sigma
	M_tau = prior$M_tau
	a_rho = prior$a_rho
	b_rho = prior$b_rho
	sigma2_beta = prior$sigma2_beta

	# Set up fixed parameters
	if (is.null(fixed)) { fixed = get_fixed() }
	stopifnot(class(fixed) == "gibbs_fixed")

	# Set up timers
	elapsed = list(beta = 0, eta = 0, sigma2 = 0, tau2 = 0, rho = 0)

	for (rep in 1:R) {
		if (rep %% report_period == 0) { logger("Starting rep %d\n", rep) }

		# Draw [eta | rest]
		st = Sys.time()
		if (!fixed$eta) {
			Omega1 = 1/sigma2 * crossprod(S, S)
			Omega2 = 1/tau2 * H
			Omega = symmpart(Omega1 + Omega2)
			a = as.matrix(y - mu)
			b = crossprod(S, a / sigma2)
			mean = solve(Omega, b)
			eta = as.numeric(rmvnorm_prec(1, mean, Omega))

			# Update phi based on eta draw
			phi = S %*% eta
		}
		elapsed$eta = elapsed$eta + as.numeric(Sys.time() - st, units = "secs")

		# Draw [beta | rest]
		st = Sys.time()
		if (!fixed$beta) {
			Omega1 = 1/sigma2 * crossprod(X, X)
			Omega2 = diag(1/sigma2_beta, d, d)
			Omega = symmpart(Omega1 + Omega2)
			a = as.matrix(y - phi)
			b = crossprod(X, a / sigma2)
			mean = solve(Omega, b)
			beta = as.numeric(rmvnorm_prec(1, mean, Omega))

			# Update mu based on beta draw
			mu = X %*% beta
		}
		elapsed$beta = elapsed$beta + as.numeric(Sys.time() - st, units = "secs")

		# Draw [sigma2 | rest]
		st = Sys.time()
		if (!fixed$sigma2) {
			a_star = n/2
			b_star = 1/2 * sum((y - mu - phi)^2)
			sigma2 = rtruncated(1, lo = 0, hi = M_sigma, pinvgamma, qinvgamma,
				a = a_star, b = b_star)
		}
		elapsed$sigma2 = elapsed$sigma2 + as.numeric(Sys.time() - st, units = "secs")

		# Draw [tau2 | rest]
		st = Sys.time()
		if (!fixed$tau2) {
			a_star = k/2
			b_star = 1/2 * as.numeric(crossprod(eta, H %*% eta))
			tau2 = rtruncated(1, lo = 0, hi = M_tau, pinvgamma, qinvgamma,
				a = a_star, b = b_star)
		}
		elapsed$tau2 = elapsed$tau2 + as.numeric(Sys.time() - st, units = "secs")

		# Draw [rho | rest]
		st = Sys.time()
		if (!fixed$rho) {
			C = as.numeric(crossprod(eta, A %*% eta)) / (2 * tau2)

			# C++ version of direct sampler
			ds_out = r_car_beta(n = 1, lambda = lambda, C = C,
				x_lower = 1e-6, x_upper = 1 - 1e-6,
				a_rho = a_rho, b_rho = b_rho,
				N = direct_control$N,
				max_rejections = direct_control$max_rejections)
			rho = ds_out$x
			rho_rejections = ds_out$rejections

			H = as.matrix(D - rho*A)
		}
		elapsed$rho = elapsed$rho + as.numeric(Sys.time() - st, units = "secs")

		if (rep > burn && rep %% thin == 0) {
			rep_keep = rep_keep + 1
			beta_hist[rep_keep,] = beta
			sigma2_hist[rep_keep] = sigma2
			tau2_hist[rep_keep] = tau2
			rho_hist[rep_keep] = rho
			rho_rejections_hist[rep_keep] = rho_rejections
			if (control$save_eta) {
				eta_hist[rep_keep,] = eta
			}
		}
	}

	ret = list(beta_hist = beta_hist, eta_hist = eta_hist, 
		sigma2_hist = sigma2_hist, tau2_hist = tau2_hist, rho_hist = rho_hist,
		R_keep = R_keep, elapsed = elapsed, R = R, burn = burn, thin = thin,
		rho_rejections_hist = rho_rejections_hist)
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

	df_tau2 = as.data.frame(cbind(
		mean(object$tau2_hist),
		sd(object$tau2_hist),
		t(quantile(object$tau2_hist, probs = pr))
	))
	rownames(df_tau2) = sprintf("tau2")

	df_rho = as.data.frame(cbind(
		mean(object$rho_hist),
		sd(object$rho_hist),
		t(quantile(object$rho_hist, probs = pr))
	))
	rownames(df_rho) = sprintf("rho")
	
	df = rbind(df_beta, df_sigma2, df_tau2, df_rho)
	colnames(df) = c("mean", "sd", "2.5%", "50%", "97.5%")
	return(df)
}

print.my_fit = function(x, ...)
{
	cat("Summary of fit for regression model with CAR random effect\n")
	print(summary(x))

	cat("----\n")
	printf("Total iterations R: %d   burn: %d   thin: %d   Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)
	printf("Rejections in rho step: %g  Total proposals: %g  Rejection rate: %0.4f%%\n",
		sum(x$rho_rejections_hist),
		sum(x$rho_rejections_hist) + x$R,
		100 * sum(x$rho_rejections_hist) / (sum(x$rho_rejections_hist) + x$R))

	cat("----\n")
	printf("Elapsed time (Seconds):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}

