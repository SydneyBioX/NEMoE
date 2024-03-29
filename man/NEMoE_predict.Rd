% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NEMoE_cv.R
\name{NEMoE_predict}
\alias{NEMoE_predict}
\title{Make predictions from a fitted "NEMoE" object.}
\usage{
NEMoE_predict(
  NEMoE_obj,
  X_new,
  Z_new = NULL,
  full = TRUE,
  level = "all",
  transform = TRUE,
  name_match = FALSE
)
}
\arguments{
\item{NEMoE_obj}{A NEMoE_obj object with fitted parameters.}

\item{X_new}{A list of transformed microbiome data.}

\item{Z_new}{A matrix of transformed nutrition data.}

\item{full}{Whether output the result in both gating network
and experts network.}

\item{level}{Which level of microbiome data to be use.}

\item{transform}{Whether use the same transformation in build NEMoE.}

\item{name_match}{Whether match the name of training data and testing data.}
}
\value{
A list of fitted result or a vector of predicted response.
}
\description{
This function predict new data using.
}
