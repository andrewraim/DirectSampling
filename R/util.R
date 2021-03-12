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

