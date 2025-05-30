% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penalizedGEE.R
\name{penalizedGEE}
\alias{penalizedGEE}
\title{Main function to solve penalized GEE}
\usage{
penalizedGEE(
  Y,
  X,
  id,
  wc.str,
  lambda,
  maxitr = 100,
  penalty = "SCAD",
  alpha_fixed = NULL
)
}
\arguments{
\item{Y}{A vector containing response values.}

\item{X}{A matrix of covariates (all numeric).}

\item{id}{A vector containing subject identifiers.}

\item{wc.str}{A character string specifying the working correlation structure. Currently allowed
structures: "independence", "exchangeable", "exchangeable_fixed", "ar1", "ar1_fixed", and
"unstructured".}

\item{lambda}{A single value of the tuning parameter.}

\item{maxitr}{Maximum number of iterations.}

\item{penalty}{The penalty type to be used, available options include "SCAD" and "MCP". 
"SCAD" refers to the Smoothly Clipped Absolute Deviation penalty and
"MCP" refers to the Minimax Concave Penalty.}

\item{alpha_fixed}{The known value of alpha which must be given when "exchangeable_fixed" or
"ar1_fixed" is used as wc.str.}
}
\value{
A list containing the estimates, asymptotic variance etc.
}
\description{
This function solves penalized GEE problem for a given longitudinal data,
 a specific working correlation structure and a single value of the tuning parameter.
}
