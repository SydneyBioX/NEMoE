% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rSample}
\alias{rSample}
\title{Cpp function: Sample from Categorical distribution}
\usage{
rSample(Prob)
}
\arguments{
\item{Prob}{A matrix with n row and p column. Each column represent a class.
The ith row is the probability distribution of ith sample.}
}
\value{
A matrix of n row and p column. The ith row draw one
sample from categorical distribution with probability in the ith row in Prob.
}
\description{
Cpp function: Sample from Categorical distribution
}
