% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{createCVList}
\alias{createCVList}
\title{Create list of parameters for cross validate NEMoE}
\usage{
createCVList(
  g1 = 10,
  lambda2_seq = c(0.005, 0.014, 0.016, 0.023, 0.025),
  shrink = 0.5,
  crit_eval = "all",
  crit_sel = "auc",
  track = FALSE,
  itmax_lambda = 3,
  itmax_cv = 20,
  itmax_fit = 50
)
}
\arguments{
\item{g1}{Numbers of parameters lambda1 for validation.}

\item{lambda2_seq}{A vector of candidates of lambda2.}

\item{shrink}{A number of shrinkage of selected lambda1. By default is 0.5.}

\item{crit_eval}{A vector of method for Evaluation metric of NEMoE object.
Can be chosen from statistics "AIC", "BIC", "ICL1", "ICL2", "eBIC",
"mAIC", "mBIC", "mICL1", "mICL2" and cross validation result "accuracy",
"D.square", "TPR", "TNR", "F1" and "auc".
If method = "all", all of the evaluation metric will be use.
By default is "all.}

\item{crit_sel}{Method use for select parameters.
Should be a member in crit_eval.}

\item{track}{Whether output the result of choosing parameters lambda1.}

\item{itmax_lambda}{Maximal iterations of generation lambda1.
By default is 3.}

\item{itmax_cv}{Maximal iterations of calculate cross validation metric.
By default is 20.}

\item{itmax_fit}{Maximal iterations of fitting NEMoE inside the evaluation.
By default is 50.}
}
\value{
A list contain parameters in fitting NEMoE.
}
\description{
This function create parameters that put into cvNEMoE function.
}
\examples{
cvparams = createCVList(g1 = 10, lambda2_seq = c(0.005, 0.015, 0.02))
}
\seealso{
\code{\link{NEMoE_buildFromList}},
\code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{fitNEMoE}},
\code{\link{calcCriterion}}, \code{\link{cvNEMoE}}
}
