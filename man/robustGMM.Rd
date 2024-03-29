% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustGMM.R
\name{robustGMM}
\alias{robustGMM}
\title{Fitting robust Gaussian Mixture Model}
\usage{
robustGMM(
  x,
  lambda = NULL,
  mu = NULL,
  sigma = NULL,
  beta = 0.2,
  tol = 1e-06,
  maxiter = 1000,
  verbose = FALSE,
  keep.data = FALSE
)
}
\arguments{
\item{x}{`numeric(n)` - A numeric vector of observations}

\item{lambda}{`numeric` - Vector of mixing proportions of each normal component}

\item{mu}{`numeric` - Vector of means of each normal component}

\item{sigma}{`numeric` - Vector of standard deviations of each normal component}

\item{beta}{`numeric(1)` - parameter corresponds to the beta-divergence induced likelihood}

\item{maxiter}{Max iteration}

\item{keep.data}{}
}
\description{
Fitting robust Gaussian Mixture Model
}
