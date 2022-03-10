###################################################################
# function to calculate AIC of NEMoE fitting result
###################################################################

.NEMoE_AIC <- function(NEMoE_obj, modified = FALSE, LL = "obs"){

  K <- NEMoE_obj@K
  LL_out <- NEMoE_obj@NEMoE_output$LL

  if(LL == "comp"){
    LogLikeli <- LL_out[nrow(LL_out),3]
  }else{
    LogLikeli <- LL_out[nrow(LL_out),1]
  }

  if(modified){

    NEMoE_df = calcdf(NEMoE_obj, output = "comp")
    df_gamma <- NEMoE_df$df_gamma
    df_beta <- NEMoE_df$df_beta

    AIC_est <- (-2*LogLikeli) + 2*(df_gamma + sum(df_beta))

  }else{
    NEMoE_df = calcdf(NEMoE_obj)
    df_gamma <- NEMoE_df$df_gamma
    df_beta <- NEMoE_df$df_beta

    AIC_est <- (-2*LogLikeli) + 2*(df_gamma + sum(df_beta))
  }

  return(-AIC_est)

}

.acc <- function(y_real, y_pred){

  res = sum(round(y_pred) == round(y_real))/length(y_real)
  return(res)
}

.dev <- function(y_real, y_pred){

  D_explain = sum(y_real*log(y_pred) +
                    (1 - y_real)*log(1 - y_pred))
  D_null = sum(y_real*log(0.5) + (1 - y_real)*log(0.5))

  return((D_null - D_explain)/D_null)
}

.TPR <- function(y_real, y_pred){

  idx <- which(y_real == 1)
  y_real_P = y_real[idx]
  y_pred_P = y_pred[idx]
  return(.acc(y_real_P, y_pred_P))
}

.TNR <- function(y_real, y_pred){

  idx <- which(y_real == 0)
  y_real_N = y_real[idx]
  y_pred_N = y_pred[idx]
  return(.acc(y_real_N, y_pred_N))
}

.F1 <- function(y_real, y_pred){

  y_pred1 = round(y_pred)
  Precision = .TPR(y_real, y_pred1)
  Recall = .TPR(y_pred1, y_real)

  return(2/(1/Precision + 1/Recall))

}

.auc <- function(y_real, y_pred){
  if(length(unique(y_real)) != 2){
    return(0.5)
  }else{
    predob<- ROCR::prediction(y_pred, y_real)
    perf.auc<- ROCR::performance(predob, measure = 'auc', x.measure = 'cutoff')
    return(perf.auc@y.values[[1]])
  }
}

#' @importFrom stats na.omit
.crit_wrapper <- function(y_real, y_pred,
                          cv_list = c("accuracy","D.square","auc")){

  f_list = list(.acc, .dev, .TPR, .TNR, .F1, .auc)
  implemented_crit = c("accuracy", "D.square", "TPR", "TNR", "F1", "auc")
  names(f_list) = implemented_crit
  result <- c()
  if("all" %in% cv_list){
    for(i in 1:length(f_list)){
      result[i] = f_list[[i]](y_real, y_pred)
    }
    names(result) = implemented_crit
  }else{
    cv_list1 = implemented_crit[na.omit(match(cv_list,implemented_crit))]
    f_sel <- f_list[cv_list1]
    for(i in 1:length(f_sel)){
      result[i] = f_sel[[i]](y_real, y_pred)
    }
    names(result) = cv_list1
  }

  return(result)
}

###################################################################
# function to calculate BIC of NEMoE fitting result
###################################################################

.NEMoE_BIC <- function(NEMoE_obj, modified = FALSE, type = "BIC", LL = "obs"){

  K <- NEMoE_obj@K
  L <- length(NEMoE_obj@Microbiome)
  LL_out <- NEMoE_obj@NEMoE_output$LL
  if(LL == "comp"){
    LogLikeli <- LL_out[nrow(LL_out),3]
  }else{
    LogLikeli <- LL_out[nrow(LL_out),1]
  }
  n <- nrow(NEMoE_obj@Nutrition)
  Z1 <- cbind(rep(1,n), NEMoE_obj@Nutrition)
  pi <- calcProb(Z1, NEMoE_obj@NEMoE_output$gamma)
  pi1 <- colSums(pi)
  q <- ncol(NEMoE_obj@Nutrition)
  p_L <- sapply(NEMoE_obj@Microbiome, ncol)

  if(modified){

    NEMoE_df = calcdf(NEMoE_obj, output = "comp")
    df_gamma <- NEMoE_df$df_gamma
    df_beta <- NEMoE_df$df_beta

    if(type == "eBIC"){
      BIC_est <- (-2 * LogLikeli) + df_gamma*log(n) + df_gamma*log(q)
      for(i in 1:L){
        for(j in 1:K){
          BIC_est <- BIC_est + df_beta[i,j]*(log(pi1[j]) + log(p_L[i]))
        }
      }
      BIC_est <- unname(BIC_est)
    }else{
      df_beta1 <- colSums(df_beta)*log(pi1)
      BIC_est <- (-2 * LogLikeli) + df_gamma*log(n) + sum(df_beta1)
    }

  }else{

    NEMoE_df = calcdf(NEMoE_obj)
    df_gamma <- NEMoE_df$df_gamma
    df_beta <- NEMoE_df$df_beta

    if(type == "eBIC"){
      BIC_est <- (-2 * LogLikeli) + df_gamma*(log(n) + log(q)) +
        sum(df_beta)*log(n) + sum(df_beta*log(p_L))
    }else{
      BIC_est <- (-2 * LogLikeli) + (df_gamma + sum(df_beta))*log(n)
    }
  }
  return(-BIC_est)
}

###################################################################
# function to calculate Cross Validation metrics of NEMoE object
###################################################################

.NEMoE_cv <- function(NEMoE_obj, crit_list = c("D.square", "accuracy", "auc"),
                      fold = 5, itmax = 1e2, verbose = F){

  L <- length(NEMoE_obj@Microbiome)
  q <- ncol(NEMoE_obj@Nutrition)
  n <- nrow(NEMoE_obj@Nutrition)
  num_fold <- ceiling(n/fold)
  idx_all <- seq_len(n)
  if("all" %in% crit_list){
    result = matrix(0, ncol = 6, nrow = fold)
    colnames(result) = c("accuracy", "D.square", "TPR", "TNR", "F1", "auc")
  }else{
    result = matrix(0, ncol = length(crit_list), nrow = fold)
    colnames(result) = crit_list
  }

  for(i in 1:fold){

    idx_test <- seq((i - 1)*num_fold + 1, min(n, i*num_fold))
    idx_train <- idx_all[-idx_test]

    Nutrition_sub = NEMoE_obj@Nutrition[idx_train,]
    Nutrition_test = NEMoE_obj@Nutrition[idx_test,]
    Microbiome_sub  = list()
    Microbiome_test = list()
    for(j in 1:L){
      Microbiome_sub[[j]] = NEMoE_obj@Microbiome[[j]][idx_train,]
      Microbiome_test[[j]] = NEMoE_obj@Microbiome[[j]][idx_test,]
    }
    Response_sub = NEMoE_obj@Response[idx_train]
    Response_test = NEMoE_obj@Response[idx_test]

    NEMoE_sub = NEMoE_obj
    NEMoE_sub@Nutrition = Nutrition_sub
    NEMoE_sub@Microbiome = Microbiome_sub
    NEMoE_sub@Response = Response_sub

    NEMoE_sub@params$verbose = verbose
    NEMoE_sub@params$itmax = itmax

    beta_init = NEMoE_sub@NEMoE_output$beta
    gamma_init = NEMoE_sub@NEMoE_output$gamma

    NEMoE_sub = .fitNEMoE(NEMoE_sub, beta_init = beta_init, gamma_init = gamma_init)
    Response_pred = NEMoE_predict(NEMoE_sub, Microbiome_test,
                                  Nutrition_test, full = F, transform = F)

    result[i,] = .crit_wrapper(Response_test, Response_pred,
                               crit_list)
  }
  return(result)
}


