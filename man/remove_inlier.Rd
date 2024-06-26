% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depurate.R
\name{remove_inlier}
\alias{remove_inlier}
\title{Remove spatial outliers}
\usage{
remove_inlier(
  x,
  y,
  ldist = 0,
  udist = 40,
  criteria = c("LM", "MP"),
  zero.policy = NULL
)
}
\arguments{
\item{x}{an \code{sf} points object}

\item{y}{\code{character} with the name of the variable to use for
depuration process}

\item{ldist}{\code{numeric} lower distance bound to identify neighbors}

\item{udist}{\code{numeric} upper distance bound to identify neighbors}

\item{criteria}{\code{character} with "LM" and/or "MP" for methods to
identify spatial outliers}

\item{zero.policy}{default NULL, use global option value;
if FALSE stop with error for any empty neighbors sets,
if TRUE permit the weights list to be formed with zero-length
weights vectors}
}
\description{
Removes spatial outliers using Local Moran's I statistic
  and moran scatterplot.
}
