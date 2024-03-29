% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{NEMoE_buildFromPhyloseq}
\alias{NEMoE_buildFromPhyloseq}
\title{Build NEMoE object from a phyloseq object.}
\usage{
NEMoE_buildFromPhyloseq(
  ps,
  Nutrition,
  Response,
  K = NULL,
  gatherParam = list(),
  filtParam = list(),
  transParam = list(),
  cvParams = list(),
  taxLevel = NULL,
  taxTab = NULL,
  standardize = TRUE,
  ...
)
}
\arguments{
\item{ps}{A phyloseq object of input microbiome data}

\item{Nutrition}{A dataframe of matrix of Nutrition matrix.}

\item{Response}{A numeric(factor) vector of Health outcome.}

\item{K}{A integer of number of latent class.}

\item{gatherParam}{A list of parameters for gather phyloseq object.
See \code{\link{psGather}}.}

\item{filtParam}{A list of parameters for filtering of Microbiome data.
can be set by prev (The number of non-zero proportion of features) and
var(varaince of features).}

\item{transParam}{A list of parameters for transformation of Microbiome data.
See \code{\link{compTransform}}.}

\item{cvParams}{A list of cross validation parameters.}

\item{taxLevel}{A character of selected name of taxonomic levels.}

\item{taxTab}{A dataframe of taxonomic table.}

\item{standardize}{Logical flag for x variable standardization.
Default is standardize=TRUE.}

\item{...}{Other parameters can pass to NEMoE_buildFromPhyloseq.
See \code{\link{createParameterList}}}
}
\value{
A NEMoE object.
}
\description{
Build NEMoE object from a phyloseq object.
}
\examples{
data(PD)
ps = PD$ps
Response <- PD$data_list$Response
Nutrition_data <- PD$data_list$Nutrition
NEMoE <- NEMoE_buildFromPhyloseq(ps = ps, Nutrition = Nutrition_data, Response = Response)
}
