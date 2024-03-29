% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{psGather}
\alias{psGather}
\title{Transformation of phyloseq data}
\usage{
psGather(
  ps,
  tax_list = c("Phylum", "Order", "Family", "Genus", "ASV"),
  prop_trans = T
)
}
\arguments{
\item{ps}{a phyloseq object.}

\item{tax_list}{Get taxonomic levels of data. By default, return the taxonomic of Phylum, Order, Family, Genus and ASV/OTU.}

\item{prop_trans}{a logic variable to determine whether the data are as proportion. By default is TRUE.}
}
\value{
A list of data of different taxonomic levels.
}
\description{
Transformation of phyloseq data
}
\examples{
data(PD)
ps1 <- psGather(PD$ps)
}
