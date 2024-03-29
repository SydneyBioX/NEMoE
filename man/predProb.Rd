% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{predProb}
\alias{predProb}
\title{Cpp function: Predict probability of mixture distribution}
\usage{
predProb(X, Z, y, W, V)
}
\arguments{
\item{X}{a data matrix of input in experts network.}

\item{Z}{a data matrix of input in gating network.}

\item{y}{a vector of response.}

\item{W}{parameters in experts network.}

\item{V}{parameters in gating network.}
}
\value{
A matrix of predicted probability with
predProb = pi * (1/(1 + exp(-X*W))),
pi = exp(Z * V_i)/sum(Z * V)
}
\description{
Cpp function: Predict probability of mixture distribution
}
