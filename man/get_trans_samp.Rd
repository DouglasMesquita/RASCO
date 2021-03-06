% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_sc.R
\name{get_trans_samp}
\alias{get_trans_samp}
\title{get_trans_samp}
\usage{
get_trans_samp(marg, fun, n = 1000, trunc = FALSE, method = "linear")
}
\arguments{
\item{marg}{a marginal from a INLA model.}

\item{fun}{a transformation function.}

\item{n}{a number of samples}

\item{trunc}{TRUE: x > 0.}

\item{method}{?inla.tmarginal method.}
}
\description{
Generating a sample from the transformed marginal
}
\keyword{internal}
