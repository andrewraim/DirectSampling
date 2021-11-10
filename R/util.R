#' Utilities
#' 
#' @param msg A string which can have format placeholders as in \code{sprintf}.
#' @param ... Arguments to \code{sprintf} to use in the format.
#' 
#' @details
#' \code{printf} mimics \code{printf} in C. \code{logger} prepends the message
#' with a time stamp, which is useful in tracking long programs.
#' 
#' @name Utilities
NULL

#' @name Utilities
#' @export
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

#' @name Utilities
#' @export
logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

# The following function helps us to use the property:
#  log(x + y) = log(x) + log(1 + y/x)
#             = log(x) + log1p(exp(log(y) - log(x))),
#  with log(x) and log(y) given as inputs, i.e. on the log-scale. When x and y
#  are of very different magnitudes, this is more stable when x is taken to be
#  the larger of the inputs. The most extreme case is when one of the inputs
#  might be -inf; in this case that input should be the second one.
#
#  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
logadd = function(logx, logy) {
	if (logx < 0 && logy < 0 && is.infinite(logx) && is.infinite(logy)) {
		# Is it possible to handle this case more naturally?
		return(-Inf)
	}
	logx + log1p(exp(logy - logx))
}

# Same as logadd, but for subtraction
#  log(x - y) = log(x) + log(1 - y/x)
#             = log(x) + log1p(-exp(log(y) - log(x)))
logsub = function(logx, logy) {
	if (logx < 0 && logy < 0 && is.infinite(logx) && is.infinite(logy)) {
		# Is it possible to handle this case more naturally?
		return(-Inf)
	}
	logx + log1p(-exp(logy - logx))
}
