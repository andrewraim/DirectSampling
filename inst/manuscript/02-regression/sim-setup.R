# Make a folder for each combination of simulation settings. Each folder
# contains a launch.R script which carries out this run of the simulation
# and leaves results in that folder.
library(DirectSampling)
source("../shared/util.R")

# Simulation settings to vary
rho_levels = c(0.01, 0.10, 0.40)
sigma_levels = c(0.25, 1, 5)

# Global simulation settings
n = 200
beta_true = c(5, -1)
N_sim = 500
mcmc_reps = 2000
mcmc_burn = 1000
report_period = 500

for (idx_rho in seq_along(rho_levels)) {
for (idx_sigma in seq_along(sigma_levels)) {
	# Create a folder for this setting of the simulation
	dd = sprintf("sim_rho%03d_sigma%03d", idx_rho, idx_sigma)
	if (dir.exists(dd)) {
		next
	}
	dir.create(dd)

	# Write launch.R file into folder dd
	ff = file(sprintf("%s/launch.R", dd), "w")
	fprintf(ff, "set.seed(1234)\n\n")
	fprintf(ff, "n = %d\n", n)
	fprintf(ff, "rho = rep(%f, n)\n", rho_levels[idx_rho])
	fprintf(ff, "sigma2_true = %f\n", sigma_levels[idx_sigma]^2)
	fprintf(ff, "beta_true = %s\n", print_vector(beta_true))
	fprintf(ff, "N_sim = %d\n", N_sim)
	fprintf(ff, "mcmc_reps = %d\n", mcmc_reps)
	fprintf(ff, "mcmc_burn = %d\n", mcmc_burn)
	fprintf(ff, "report_period = %d\n\n", report_period)
	fprintf(ff, "source(\"../sim-rep.R\", chdir = TRUE)\n")
	fprintf(ff, "save.image(\"results.Rdata\")\n")
	close(ff)
}
}

