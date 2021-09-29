###################################################################
# function to predict output of new data using fitted NEMoE object
###################################################################
#' Make predictions from a fitted "NEMoE" object.
#' @description This function predict new data using.
#' @param NEMoE A NEMoE object with fitted parameters.
#' @param X_new A list of transformed microbiome data.
#' @param Z_new A matrix of transformed nutrition data.
#' @param full Whether output the result in both gating network
#'  and experts network.
#' @param level Which level of microbiome data to be use.
#' @param transform Whether use the same transformation in build NEMoE.
#' @param name_match Whether match the name of training data and testing data.
#' @return A list of fitted result or a vector of predicted response.
#' @useDynLib NEMoE
#' @export

NEMoE_predict <- function(NEMoE, X_new, Z_new = NULL,
                          full = TRUE, level = "all",
                          transform = TRUE,
                          name_match = FALSE){

  K <- NEMoE@K
  beta <- NEMoE@NEMoE_output$beta
  gamma <- NEMoE@NEMoE_output$gamma
  L <- length(beta)
  n <- nrow(X_new[[1]])
  .transformation = NEMoE@.transformation
  method = .transformation$method

  if(length(.transformation$mu_Z) && length(.transformation$sd_Z)){
    scale_Z = T
  }else{
    scale_Z = F
  }

  if(is.null(Z_new)){
    probs = rep(1/K, K)
  }else{
    if(name_match){
      Z_new <- .mathcNEMoE(NEMoE@Nutrition, Z_new)
    }
    if(scale_Z){
      mu_Z <- .transformation$mu_Z
      sd_Z <- .transformation$mu_Z
      Z_new <- (Z_new - mu_Z)/sd_Z
    }
    Z1 = cbind(rep(1,n), Z_new)
    probs = calcProb(Z1, gamma)
  }

  if(length(.transformation$mu_X) && length(.transformation$sd_X)){
    scale_X = T
  }else{
    scale_X = F
  }

  if(level == "all"){
    probs_sub <- list()
    y = matrix(0, nrow = n, ncol = L)
    for(i in 1:L){
      X_temp <- X_new[[i]]
      if(name_match){
        X_temp <- .mathcNEMoE(NEMoE@Microbiome[[i]], X_temp)
      }
      if(transform){
        X_temp <- compTransform(X_temp, method = method, scale = F)
      }
      if(scale_X){
        mu_X_temp <- .transformation$mu_X[[i]]
        sd_X_temp <- .transformation$sd_X[[i]]
        X_temp <- (X_temp - mu_X_temp)/sd_X_temp
      }
      X_temp1 <- cbind(rep(1,n),X_temp)
      probs_temp = X_temp1 %*% beta[[i]]
      probs_temp = matrix(sapply(probs_temp, .logistic), nrow = nrow(probs_temp))
      probs_sub[[i]] <- probs_temp
      y[,i] = rowSums(probs_temp*probs)
    }
    names(probs_sub) = names(NEMoE@taxLevel)

  }else{
    y = matrix(0, nrow = n, ncol = 1)
    beta_i = beta[[level]]
    X_new <- .mathcNEMoE(NEMoE@Microbiome[[level]], X_new)
    if(scale_X){
      mu_X_temp <- .transformation$mu_X[[level]]
      sd_X_temp <- .transformation$sd_X[[level]]
      X_new1 <- (X_new - mu_X_temp)/sd_X_temp
    }
    X_temp1 <- cbind(rep(1,n), X_new1)
    probs_temp <- X_temp1 %*% beta_i
    probs_temp = matrix(sapply(probs_temp, .logistic), nrow = nrow(probs_temp))
    probs_sub = list()
    probs_sub[[level]] = probs_temp
    y = rowSums(probs_temp * probs)
  }

  if(full){
    return(list(gating_prob = probs, experts_prob = probs_sub, output = y))
  }else{
    return(rowMeans(y))
  }
}

.lambda1Generator <- function(NEMoE, lambda2, itmax = 5){

  L = length(NEMoE@Microbiome)
  y <- NEMoE@Response
  NEMoE_temp = NEMoE
  NEMoE_temp@params$lambda2 = lambda2
  gamma_init = NEMoE@NEMoE_output$gamma
  NEMoE_temp@params$verbose = F
  lambda1_M = c()
  for(i in 1:L){
    lambda1_max = max(abs(stats::cor(NEMoE@Microbiome[[i]], y)))
    repeat{
      NEMoE_temp@Microbiome = list(NEMoE@Microbiome[[i]])
      NEMoE_temp@params$lambda1 = lambda1_max
      NEMoE_temp@params$itmax = itmax
      NEMoE_temp = .fitNEMoE(NEMoE_temp, beta_init = NULL, gamma_init = gamma_init)
      df_beta = .calcdf(NEMoE_temp@NEMoE_output$beta[[1]])
      if(!df_beta){
        lambda1_max = lambda1_max*0.8
      }else{
        break()
      }
    }
    lambda1_M[i] = lambda1_max
  }
  return(lambda1_M)

}

###################################################################
# function to calculate evaluation metric of NEMoE object
###################################################################
#' Evaluate fitted NEMoE obejct
#' @description This function calculate evaluation of NEMoE object.
#' @param NEMoE a NEMoE of object with fitted parameters.
#' @param crit_list a list of evaluation metric.
#' @param ... other parameters can be put into calcCriterion.
#' @return A vector of evaluation result.
#' @examples
#' data(NEMoE_example)
#' calcCriterion(NEMoE_example)
#' @export
#' @seealso \code{\link{createCVList}}


calcCriterion <- function(NEMoE, crit_list= c("AIC", "BIC", "ICL1"),
                          ...){

  f_AIC <- .NEMoE_AIC
  f_BIC <- .NEMoE_BIC
  f_ICL1 <- purrr::partial(.NEMoE_AIC, LL = "comp")
  f_ICL2 <- purrr::partial(.NEMoE_BIC, LL = "comp")
  f_eBIC <- purrr::partial(.NEMoE_BIC, type = "eBIC")
  f_mAIC <- purrr::partial(.NEMoE_AIC, modified = TRUE)
  f_mBIC <- purrr::partial(.NEMoE_BIC, modified = TRUE)
  f_mICL1 <- purrr::partial(.NEMoE_AIC, modified = TRUE, LL = "comp")
  f_mICL2 <- purrr::partial(.NEMoE_BIC, modified = TRUE, LL = "comp")

  if("all" %in% crit_list){
    crit_list <- c("AIC", "BIC", "ICL1", "ICL2", "eBIC", "mAIC", "mBIC",
                   "mICL1", "mICL2", "accuracy", "D.square", "TPR",
                   "TNR", "F1", "auc")
  }

  f_list = list(f_AIC, f_BIC, f_ICL1, f_ICL2, f_eBIC, f_mAIC, f_mBIC,
                f_mICL1, f_mICL2)
  stat_list <- c("AIC", "BIC", "ICL1", "ICL2", "eBIC", "mAIC", "mBIC",
                 "mICL1", "mICL2")
  names(f_list) <- stat_list
  cv_list <- c("accuracy", "D.square", "TPR", "TNR", "F1", "auc")

  stat_sel <- stat_list[na.omit(match(crit_list, stat_list))]
  f_stat <- f_list[stat_sel]
  cv_sel <- cv_list[na.omit(match(crit_list, cv_list))]

  result_stat <- c()
  if(length(stat_sel)){
    for(i in 1:length(stat_sel)){
      result_stat[i] <- f_stat[[i]](NEMoE)
    }
  }
  names(result_stat) <- stat_sel

  result_cv_m <- c()
  if(length(cv_sel)){
    result_cv <- .NEMoE_cv(NEMoE, crit_list = cv_sel, ...)
    result_cv_m <- colMeans(result_cv)
  }

  result <- matrix(c(result_stat, result_cv_m), nrow = 1)
  colnames(result) <- c(stat_sel, cv_sel)

  return(result)
}

.cvNEMoE <- function(NEMoE, g1 = 10, itmax_lambda = 3,
                    itmax_fit = 50, itmax_cv = 20,
                    crit_list = c("all"),
                    lambda1_M = NULL, ...){

  L <- length(NEMoE@Microbiome)
  p_L <- sapply(NEMoE@Microbiome, ncol)
  n <- nrow(NEMoE@Nutrition)
  K <- NEMoE@K

  lambda2 = NEMoE@params$lambda2

  if(is.null(lambda1_M)){
    lambda1_M <- .lambda1Generator(NEMoE, lambda2, itmax_lambda)
  }

  lambda1_m <- log(p_L)/(K*n)

  gamma_init <- NEMoE@NEMoE_output$gamma
  result = list()

  if("all" %in% crit_list){
    crit_list <- c("AIC", "BIC", "ICL1", "ICL2", "eBIC", "mAIC", "mBIC",
                   "mICL1", "mICL2", "accuracy", "D.square", "TPR",
                   "TNR", "F1", "auc")
  }

  for(i in 1:L){
    param_temp <- matrix(0, nrow = g1, ncol = 2)
    result_temp <- matrix(0, nrow = g1, ncol = length(crit_list))
    C <- (log(lambda1_M[i]) - log(lambda1_m[i])) / (g1 - 1)
    lambda1_seq <- lambda1_m[i] * exp((seq(1,g1) - 1)*C)
    NEMoE_temp <- NEMoE
    NEMoE_temp@Microbiome <- list(NEMoE_temp@Microbiome[[i]])
    NEMoE_temp@params$verbose = FALSE
    NEMoE_temp@params$itmax = itmax_fit
    NEMoE_temp@params$early_stop = TRUE

    for(j in 1:g1){
      message("Evaluate on lambda1 = ", lambda1_seq[j],"...")
      NEMoE_temp@params$lambda1 = lambda1_seq[j]
      NEMoE_temp <- .fitNEMoE(NEMoE_temp,
                              beta_init = NULL,
                              gamma_init = gamma_init)
      param1 = paste(NEMoE_temp@params$lambda1, collapse = ",")
      param2 = paste(NEMoE_temp@params$lambda2, collapse = ",")
      param_temp[j,] <- c(param1, param2)
      res_temp <- calcCriterion(NEMoE_temp, crit_list = crit_list,
                                itmax = itmax_cv, ...)
      result_temp[j,] <- res_temp
    }
    param_df <- as.data.frame(param_temp)
    colnames(param_df) <- c("lambda1", "lambda2")
    result_df <- as.data.frame(result_temp)
    colnames(result_df) <- colnames(res_temp)
    result[[i]] <- cbind(param_df, result_df)
  }
  return(result)
}

.chooselambda1 <- function(result, crit_sel = "auc", shrink = 0.5){

  lambda1_sel <- c()
  idx_sel <- which.max(result[,crit_sel])
  lambda1_sel <- as.numeric(result[idx_sel,1]) * shrink

  return(lambda1_sel)
}

###########################################################################
# function to Cross validate NEMoE with different tunning parameters lambda
##########################################################################
#' Parameters tunning in NEMoE
#' @description This function calculate evaluation of NEMoE object.
#' @param NEMoE a NEMoE of object without cross validated parameters.
#' @param verbose  A logical input indicating whether the intermediate
#' steps will be printed.
#' @param ... Other parameters that can be passed to cvNEMoE.
#' @return A NEMoE object with different lambdas.
#' @examples
#' data(NEMoE_example)
#' NEMoE_example = cvNEMoE(NEMoE_example, verbose = F)
#' @export
#' @seealso \code{\link{createCVList}}

cvNEMoE <- function(NEMoE, verbose = T, ...){

  lambda2_seq <- NEMoE@cvParams$lambda2
  g1 <- NEMoE@cvParams$g1
  shrink <- NEMoE@cvParams$shrink
  itmax_cv <- NEMoE@cvParams$itmax_cv
  itmax_fit <- NEMoE@cvParams$itmax_fit
  itmax_lambda <- NEMoE@cvParams$itmax_lambda
  crit_eval <-  NEMoE@cvParams$crit_eval
  crit_sel <- NEMoE@cvParams$crit_sel
  track <- NEMoE@cvParams$track
  gamma_init <- NEMoE@NEMoE_output$gamma

  if("all" %in% crit_eval){
    crit_eval <- c("AIC", "BIC", "ICL1", "ICL2", "eBIC", "mAIC", "mBIC",
                   "mICL1", "mICL2", "accuracy", "D.square", "TPR",
                   "TNR", "F1", "auc")
  }

  result_lambda1 <- list()
  lambda1_sel <- list()
  result_lambda2 <- matrix(0, nrow = length(lambda2_seq),
                           ncol = length(crit_eval))
  param_temp <- matrix(0, nrow = length(lambda2_seq), ncol = 2)

  for(i in 1:length(lambda2_seq)){
    message("Evaluate on lambda2 = ", lambda2_seq[i], ".....")
    NEMoE_temp <- NEMoE
    NEMoE_temp@params$lambda2 = lambda2_seq[i]
    result_lambda1[[i]] <- .cvNEMoE(NEMoE_temp, g1, itmax_lambda,
                                    itmax_fit, itmax_cv,
                                    crit_eval)

    lambda1_sel[[i]] <- sapply(result_lambda1[[i]],
                               .chooselambda1,
                               crit_sel = crit_sel,
                               shrink = shrink)

    NEMoE_temp@params$lambda1 <- lambda1_sel[[i]]
    NEMoE_temp@params$verbose = FALSE
    NEMoE_temp@params$itmax = itmax_fit
    NEMoE_temp@params$early_stop = TRUE

    NEMoE_temp <- .fitNEMoE(NEMoE_temp,
                            beta_init = NULL,
                            gamma_init = gamma_init)

    result_temp <- calcCriterion(NEMoE = NEMoE_temp, crit_list = crit_eval,
                                 itmax = itmax_cv, ...)

    result_lambda2[i,] <- result_temp

    param1 = paste(NEMoE_temp@params$lambda1, collapse = ",")
    param2 = paste(NEMoE_temp@params$lambda2, collapse = ",")
    param_temp[i,] <- c(param1, param2)
  }
  param_df <- as.data.frame(param_temp)
  colnames(param_df) <- c("lambda1", "lambda2")
  result_df <- as.data.frame(result_lambda2)
  colnames(result_df) <- colnames(result_temp)
  result_lambda2 <- cbind(param_df, result_df)

  idx_sel <- which.max(result_lambda2[,crit_sel])
  lambda1_choose <- lambda1_sel[[idx_sel]]
  lambda2_choose <- lambda2_seq[idx_sel]

  if(track){

    cv_result <- list(result_lambda1 = result_lambda1,
                      result_lambda2 = result_lambda2,
                      lambda1_choose = lambda1_choose,
                      lambda2_choose = lambda2_choose)
  }else{
    cv_result <- list(result_lambda2 = result_lambda2,
                      lambda1_choose = lambda1_choose,
                      lambda2_choose = lambda2_choose)
  }
  NEMoE@cvResult <- cv_result

  return(NEMoE)
}
