% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{confint.medoutcon}
\alias{confint.medoutcon}
\title{Confidence intervals for natural/interventional (in)direct effect estimates (This function is directly from the medoutcon package)}
\usage{
\method{confint}{medoutcon}(object, parm = seq_len(object$theta), level = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{medoutcon}, as produced by invoking
\code{\link{medoutcon}}, for which a confidence interval is to be computed.}

\item{parm}{A \code{numeric} vector indicating indices of \code{object$est}
for which to return confidence intervals.}

\item{level}{A \code{numeric} indicating the level of the confidence
interval to be computed.}

\item{...}{Other arguments. Not currently used.}
}
\description{
Compute confidence intervals for objects of class \code{medoutcon}, which
contain estimates produced by \code{\link{medoutcon}}.
}
\keyword{internal}
