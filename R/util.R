#' Utilities
#' 
#' @param msg A string
#' 
#' @details
#' \code{printf} mimics \code{printf} in C. \code{logger} prepends the message
#' with a time stamp, which is useful in tracking long programs.
#' 
#' @name Util
NULL

#' @name Util
#' @export
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

#' @name Util
#' @export
logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

