% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NEMoE_example}
\alias{NEMoE_example}
\title{An example of fitted NEMoE object}
\format{
An object of class \code{NEMoE} of length 1.
}
\usage{
data(NEMoE_example, package = 'NEMoE')
}
\description{
An example of fitted NEMoE object
Generate from following code
data(PD)
Microbiome = PD$data_list$Microbiome
Nutrition = PD$data_list$Nutrition
Response = PD$data_list$Response
NEMoE = NEMoE_buildFromList(Microbiome, Nutrition, Response,
lambda1 = c(0.005, 0.014, 0.016, 0.023, 0.025),lambda2 = 0.02,
alpha1 = 0.5, alpha2 = 1, adapt = T, btr = T, stop_all = T, itmax = 1e3,
verbose = T, beta_max = 1e2, early_stop = F)
NEMoE = fitNEMoE(NEMoE, num_restart = 0)
}
\keyword{datasets}
