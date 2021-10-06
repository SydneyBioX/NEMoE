###################################################################
# function to build NEMoE object from phyloseq object
###################################################################
#' Build NEMoE object from a phyloseq object.
#' @param ps A phyloseq object of input microbiome data
#' @param Nutrition A dataframe of matrix of Nutrition matrix.
#' @param Response A numeric(factor) vector of Health outcome.
#' @param K A integer of number of latent class.
#' @param gatherParam A list of parameters for gather phyloseq object.
#' See \code{\link{psGather}}.
#' @param filtParam A list of parameters for filtering of Microbiome data.
#' can be set by prev (The number of non-zero proportion of features) and
#' var(varaince of features).
#' @param transParam A list of parameters for transformation of Microbiome data.
#' See \code{\link{compTransform}}.
#' @param cvParams A list of cross validation parameters.
#' @param taxLevel A character of selected name of taxonomic levels.
#' @param taxTab A dataframe of taxonomic table.
#' @param TSS A logical variable indicate whether use TSS normalization.
#' @param ... Other parameters can pass to NEMoE_buildFromPhyloseq.
#' See \code{\link{createParameterList}}
#' @return A NEMoE object.
#' @examples
#' data(PD)
#' ps = PD$ps
#' Response <- PD$data_list$Response
#' Nutrition_data <- PD$data_list$Nutrition
#' NEMoE <- NEMoE_buildFromPhyloseq(ps = ps, Nutrition = Nutrition_data, Response = Response)
#' @export
NEMoE_buildFromPhyloseq <- function(ps, Nutrition, Response, K = NULL,
                                    gatherParam = list(), filtParam = list(),
                                    transParam = list(), cvParams = list(),
                                    taxLevel = NULL, taxTab = NULL, ...,
                                    TSS = T){
  if(is.null(taxLevel)){
    taxLevel = c("Phylum","Order","Family","Genus","ASV")
  }
  gatherParam$tax_list = taxLevel

  Microbiome <- do.call(psGather, c(ps, gatherParam))
  L <- length(Microbiome)
#  for(i in 1:L){
#    Microbiome[[i]] <- matrix(Microbiome[[i]], nrow = nrow(Microbiome[[i]]))
#  }
  .transformation = list(method = "none", mu_X = list(), sd_X = list())

  if(!length(filtParam)){
    filtParam = list(prev = 0.7, var = 5e-5)
  }

  if(!length(transParam)){
    transParam = list(method = "comp", scale = T)
  }

  ps <- phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})

  for(i in 1:L){

    if(taxLevel[i] %in% c("ASV","OTU", "asv", "otu")){
      id_out = TRUE
      id0 = 1:ncol(Microbiome[[i]])
    }else{
      id_out = FALSE
    }

    if(!is.null(filtParam$prev)){
      q_filt = purrr::partial(stats::quantile, probs = filtParam$prev)
      M_temp <- .filterComp(Microbiome[[i]], q_filt, 0, id_out)
      Microbiome[[i]] <- M_temp$X
      if(id_out){
        id0 <- id0[M_temp$id]
      }
    }

    if(!is.null(filtParam$var)){
      M_temp <- .filterComp(Microbiome[[i]], stats::var, filtParam$var, id_out)
      Microbiome[[i]] <- M_temp$X
      if(id_out){
        id0 <- id0[M_temp$id]
      }
    }
  }

  for(i in 1:L){
    Microbiome[[i]] <- do.call(compTransform,c(list(Microbiome[[i]]), transParam))
    .transformation$method = transParam$method
    .transformation$mu_X[[i]] = attr(Microbiome[[i]], "scaled:center")
    .transformation$sd_X[[i]] = attr(Microbiome[[i]], "scaled:scale")
  }

  if(is.data.frame(Nutrition) || is.matrix(Nutrition)){
    Nutrition = as.matrix(Nutrition)
  }

  .transformation$mu_Z = attr(Nutrition, "scaled:center")
  .transformation$sd_Z = attr(Nutrition, "scaled:scale")

  taxTab = as.data.frame(ps@tax_table)[id0,]


  if(is.factor(Response)){
    ResponseLevel = levels(Response)
    Response = as.numeric(Response) - 1
  }else{
    ResponseLevel = NULL
    Response = as.numeric(Response)
  }

  params = createParameterList(...)

  if(!length(cvParams)){
    cvParams = createCVList()
  }

  if(is.null(K)){
    K = 2
  }

  NEMoE = .NEMoE(Microbiome = Microbiome, Nutrition = Nutrition,
                 Response = Response, params = params,
                 cvParams = cvParams, K = K, taxLevel = taxLevel,
                 taxTab = taxTab, ResponseLevel = ResponseLevel,
                 .transformation = .transformation)

  return(NEMoE)

}

###################################################################
# function to build NEMoE object from List, dataframe or matrix
###################################################################
#' Build NEMoE object from list, dataframe or matrix
#' @description This function build NEMoE object from List, dataframe or matrix.
#' @param Microbiome A list of transformed microbiome
#' data of different taxonomic level. If it is a dataframe or matrix, will
#' convert to a list automatically.
#' @param Nutrition A dataframe of matrix of Nutrition matrix.
#' @param Response A numeric(factor) vector of Health outcome.
#' @param K A integer of number of latent class.
#' @param cvParams A list of cross validation parameters.
#' @param taxLevel A character of selected name of taxonomic levels.
#' @param taxTab A dataframe of taxonomic table.
#' @param ... Other parameters can pass to NEMoE_buildFromList.
#' See \code{\link{createParameterList}}
#' @return A NEMoE object.
#' @examples
#' data(PD)
#' data_list = PD$data_list
#' Microbiome = data_list$Microbiom
#' Nutrition = data_list$Nutrition
#' Response = data_list$Response
#' NEMoE <- NEMoE_buildFromList(Microbiome = Microbiome, Nutrition = Nutrition, Response = Response)
#' @export
#'
NEMoE_buildFromList <- function(Microbiome, Nutrition, Response,
                                K = NULL, cvParams = list(),
                                taxLevel = NULL, taxTab = NULL, ...){

  .transformation = list(method = "none", mu_X = list(), sd_X = list())

  if(is.data.frame(Microbiome) || is.matrix(Microbiome)){
    Microbiome = list(as.matrix(Microbiome))
  }
  L <- length(Microbiome)

  if(is.data.frame(Nutrition) || is.matrix(Nutrition)){
    Nutrition = as.matrix(Nutrition)
  }

  if(is.factor(Response)){
    ResponseLevel = levels(Response)
    Response = as.numeric(Response) - 1
  }else{
    ResponseLevel = NULL
    Response = as.numeric(Response)
  }

  params = createParameterList(...)

  if(!length(cvParams)){
    cvParams = createCVList()
  }

  if(is.null(K)){
    K = 2
  }

  NEMoE = .NEMoE(Microbiome = Microbiome, Nutrition = Nutrition,
                 Response = Response, params = params,
                 cvParams = cvParams, K = K, taxLevel = taxLevel,
                 taxTab = taxTab, ResponseLevel = ResponseLevel,
                 .transformation = .transformation)
  return(NEMoE)

}

###################################################################
# function to create list of parameters for fitting NEMoE
###################################################################
#' Create list of parameters for fitting NEMoE
#' @description This function create parameters that put into fitNEMoE function.
#' @param lambda1 Penalty regularizer for the experts.
#' @param lambda2 Penalty regularizer for the gating network.
#' @param alpha1 Elastic net penalty value for experts.
#' @param alpha2 Elastic net penalty value for gating network.
#' @param stop_all Method of stop criterion. If stop_all = TRUE means that either
#'  coefficient or loss function converge. If stop_all = FALSE means that both
#'  coefficient and loss function converge.
#' @param EM_alg Method for Expecation maximization update.
#'  Can be chosen from "EM", "CEM", "GEM", "SEM", "SAEM", "GEM".
#'  By default is "EM".
#' @param itmax Maximium number of iteration in fitting NEMoE. By default is 100.
#' @param itmin Minimium number of iteration in fitting NEMoE. By default is 3.
#' @param btr Whether use backtracking during the iteration. By default is TRUE.
#' @param init Method for initialization.
#' Can be chosen from "rand", "kmeans" and "glmnet".
#' If init="rand" will use a dirichlet distribution initialing the latent class.
#' If init="kmeans" the latent class will initialized using kmeans
#' clustering of input for gating network.
#' If init="glmnet" will use a equal probability
#' with lasso as its corresponding coefficients in experts network.
#' @param beta_max Maximal of coefficient. By default is 10.
#' @param adapt whether use adaptive mode in optimization. By default is TRUE.
#' @param verbose A logical input indicating whether the intermediate
#' steps will be printed.
#' @param early_stop A logical input indicate whether early stop when one of
#' the fitted latent class have zero variables selected (to save time).
#' By default is TRUE.
#' @return A list contain parameters in fitting NEMoE.
#' @examples
#' params = createParameterList(lambda1 = c(0.005, 0.01, 0.02, 0.025))
#' @export
#' @seealso \code{\link{NEMoE_buildFromList}},
#' \code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{fitNEMoE}}

createParameterList <- function(lambda1 = 0.02, lambda2 = 0.015,
                                alpha1 = 0.5, alpha2 = 0.5,
                                beta_max = 10, EM_alg = "EM", init = "kmeans",
                                itmax = 1e2, itmin = 1, adapt = TRUE,
                                btr = TRUE, stop_all = TRUE, verbose = TRUE,
                                early_stop = FALSE){

  params = list(lambda1 = lambda1, lambda2 = lambda2,
                alpha1 = alpha1, alpha2 = alpha2,
                beta_max = beta_max, EM_alg = EM_alg, init = init,
                itmax = itmax, itmin = itmin, adapt = adapt,
                btr = btr, stop_all = stop_all, verbose = verbose,
                early_stop = early_stop)

  return(params)
}

###################################################################
# function to create list of parameters for cross validate NEMoE
###################################################################
#' Create list of parameters for cross validate NEMoE
#' @description This function create parameters that put into cvNEMoE function.
#' @param g1 Numbers of parameters lambda1 for validation.
#' @param lambda2_seq A vector of candidates of lambda2.
#' @param shrink A number of shrinkage of selected lambda1. By default is 0.5.
#' @param crit_eval A vector of method for Evaluation metric of NEMoE object.
#'  Can be chosen from statistics "AIC", "BIC", "ICL1", "ICL2", "eBIC",
#'  "mAIC", "mBIC", "mICL1", "mICL2" and cross validation result "accuracy",
#'  "D.square", "TPR", "TNR", "F1" and "auc".
#'  If method = "all", all of the evaluation metric will be use.
#'  By default is "all.
#' @param crit_sel Method use for select parameters.
#' Should be a member in crit_eval.
#' @param track Whether output the result of choosing parameters lambda1.
#' @param itmax_lambda Maximal iterations of generation lambda1.
#' By default is 3.
#' @param itmax_cv Maximal iterations of calculate cross validation metric.
#' By default is 20.
#' @param itmax_fit Maximal iterations of fitting NEMoE inside the evaluation.
#' By default is 50.
#' @return A list contain parameters in fitting NEMoE.
#' @examples
#' cvparams = createCVList(g1 = 10, lambda2_seq = c(0.005, 0.015, 0.02))
#' @export
#' @seealso \code{\link{NEMoE_buildFromList}},
#' \code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{fitNEMoE}},
#' \code{\link{calcCriterion}}, \code{\link{cvNEMoE}}

createCVList <- function(g1 = 10,
                         lambda2_seq = c(0.005, 0.014, 0.016, 0.023, 0.025),
                         shrink = 0.5, crit_eval = "all", crit_sel = "auc",
                         track = FALSE, itmax_lambda = 3,
                         itmax_cv = 20, itmax_fit = 50){

  cvParams = list(g1 = g1, lambda2_seq = lambda2_seq, shrink = shrink,
                  crit_eval = crit_eval, crit_sel = crit_sel,
                  track = track, itmax_lambda = itmax_lambda,
                  itmax_cv = itmax_cv, itmax_fit = itmax_fit)
  return(cvParams)
}



###################################################################
# function to set parameters list of NEMoE object
###################################################################
#' Set parameters of fitting NEMoE object
#' @description This function set parameters in NEMoE object.
#' @param NEMoE A NEMoE object.
#' @param ... Other parameters can pass to setParam.
#'
#' @return A NEMoE with user input parameters.
#' @examples
#' data(NEMoE_example)
#' NEMoE_example <- setParam(NEMoE_example, lambda1 = 0.03)
#' @export
#' @seealso \code{\link{createParameterList}}

setParam <- function(NEMoE, ...){

  params <- createParameterList(...)
  NEMoE@params <- params

  return(NEMoE)

}

###################################################################
# function to get fitted result of NEMoE object
###################################################################
#' Get coefficients of fitted NEMoE object
#' @description This function get coefficients in NEMoE object.
#' @param NEMoE A NEMoE object.
#'
#' @return A list contain fitted result.
#' @export
#' @seealso \code{\link{fitNEMoE}}

getCoef <- function(NEMoE){

  gating <- NEMoE@NEMoE_output$gamma
  experts <- NEMoE@NEMoE_output$beta

  return(list(coef.gating = gating, coef.experts = experts))

}


###################################################################
# function to get fitted Likelihood of NEMoE object
###################################################################
#' Get Likelihood of fitted NEMoE object
#' @description This function get Likelihood result in NEMoE object.
#' @param NEMoE A NEMoE object.
#'
#' @return A dataframe contain fitted result.
#' @export
#' @seealso \code{\link{fitNEMoE}}

getLL <- function(NEMoE){

  return(NEMoE@NEMoE_output$LL)

}
