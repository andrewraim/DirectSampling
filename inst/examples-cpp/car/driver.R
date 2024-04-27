library(CARBayes)
library(CARBayesdata)
library(sp)
library(spdep)
# library(spatialreg)
library(tidyverse)
library(Matrix)

source("gibbs.R")

stopifnot(packageVersion("CARBayes") == "1.6")
stopifnot(packageVersion("CARBayesdata") == "2.2")

# ----- Settings for MCMC chain ----
mcmc_draws = 5000
mcmc_burn = 2000
mcmc_thin = 10
mcmc_report = 100

# ----- Prepare example data-----
data("pricedata", package = "CARBayesdata")
data("GGHB.IG", package = "CARBayesdata")

head(pricedata)
pricedata = pricedata %>% mutate(logprice = log(pricedata$price))
pricedata_sp = merge(x = GGHB.IG, y = pricedata, by = "IG", all.x = FALSE)
pricedata_sp = spTransform(pricedata_sp , CRS("+proj=longlat +datum=WGS84 +no_defs"))

W_nb = poly2nb(pricedata_sp, row.names = pricedata_sp@data$IG)
W_list = nb2listw(W_nb, style = "B")
A = as(W_list, "CsparseMatrix")

y = pricedata_sp@data$logprice
X = model.matrix(~ log(crime) + rooms + sales + factor(type) + log(driveshop),
	data = pricedata_sp@data)
S = Diagonal(n = ncol(A))

# ----- Gibbs sampler with C++ direct sampling step -----
init = get_init(d = ncol(X), k = ncol(S), rho = 0.5)
prior = get_prior(M_sigma = 1000, M_tau = 1000, a_rho = 1, b_rho = 1,
	sigma2_beta = 1000, d = ncol(X))
direct_control = get_direct_control(N = 30, max_rejections = 5000)
control = get_gibbs_control(save_eta = FALSE, R = mcmc_draws, burn = mcmc_burn,
	thin = mcmc_thin, report_period = mcmc_report, direct = direct_control)
fixed = get_fixed()

set.seed(1234)
gibbs_out = gibbs_sampler(y, X, S, A, init, prior, control, fixed = fixed)
print(gibbs_out)
