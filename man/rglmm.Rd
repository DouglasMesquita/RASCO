% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_glmm.R
\name{rglmm}
\alias{rglmm}
\title{rglmm}
\usage{
rglmm(
  beta = c(0.1, -0.1),
  tau = 1,
  family = "gaussian",
  confounding = "none",
  neigh,
  scale = TRUE,
  control_family = NULL
)
}
\arguments{
\item{beta}{coefficients parameters.}

\item{tau}{precision for ICAR.}

\item{family}{"gaussian", "binomial" or "poisson".}

\item{confounding}{"none", "linear", "quadratic" or "cubic".}

\item{neigh}{neighborhood structure. A \code{SpatialPolygonsDataFrame} object.}

\item{scale}{scale covariates? TRUE or FALSE. See ?scale for more information.}

\item{control_family}{a list with, at least, a \code{invlink}.}
}
\description{
Generating data from Generalized Linear Mixed models
}
