# Mimic printf function in C
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

# Mimic fprintf function in C
fprintf = function(file, msg, ...)
{
	cat(sprintf(msg, ...), file = file)
}

# printf with timestamp
logger = function(fmt, ..., dt_fmt = "%Y-%m-%d %H:%M:%S", join = " - ")
{
    sys.time = format(Sys.time(), format = dt_fmt)
    cat(sys.time, join, sprintf(fmt, ...), sep = "")
}

# This is from the LearnBayes package
rtruncated = function (n, lo, hi, pf, qf, ...) {
	qf(pf(lo, ...) + runif(n) * (pf(hi, ...) - pf(lo, ...)), ...)
}
