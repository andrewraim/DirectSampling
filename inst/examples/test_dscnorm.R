library(DirectSampling)

n = 200000
sigma = 0.5
y = r_dscnorm(n, sigma)

tab = table(y)
vals = as.integer(names(tab))
freqs = as.integer(tab)

x_seq = seq(min(y), max(y))
f_seq_unnormalized = d_dscnorm_unnormalized(x_seq, sigma)
f_seq = f_seq_unnormalized / sum(f_seq_unnormalized)

my_xlim = range(y)
my_ylim = c(0, max(f_seq, freqs/n))

# Compare draws of discrete normal to its density (to check our code)
my_xlim = range(x_seq)
my_ylim = c(0, max(f_seq, freqs/n))
plot(vals, freqs / n, xlim = my_xlim, ylim = my_ylim)
points(x_seq, f_seq, pch = 3)

# Compare draws of discrete normal to continuous normal
my_xlim = range(x_seq)
my_ylim = c(0, max(f_seq, dnorm(0, 0, sigma)))
plot(vals, freqs / n, xlim = my_xlim, ylim = my_ylim)
curve(dnorm(x, 0, sigma), add = TRUE)


# Let's see what happens if we just truncate
sigma = 0.5
x_seq = seq(-10, 10)
f_seq_unnormalized = d_dscnorm_unnormalized(x_seq, sigma)
cbind(x_seq, f_seq_unnormalized)

# Normalize to see how insignificant the points away from zero really are
f_seq = f_seq_unnormalized / sum(f_seq_unnormalized)
cbind(x_seq, f_seq)
