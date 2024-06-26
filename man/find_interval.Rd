% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{find_interval}
\alias{find_interval}
\title{Find Interval}
\usage{
find_interval(x, cutpoints)
}
\arguments{
\item{x}{A number}

\item{cutpoints}{A sorted vector}
}
\description{
Find Interval
}
\details{
Elements of cutpoints represent endpoints
of adjacent intervals \eqn{[c_0,c_1)}, \eqn{[c_1,c_2)} ...,
\eqn{[c_{k},c_{k+1})}. Return the index \code{i} such that \code{x} is in 
\eqn{[c_i,c_{i+1})}; or return \code{-1} if \eqn{x < c_0}
or \code{k+1} if \eqn{x > c_{k+1}}
}
