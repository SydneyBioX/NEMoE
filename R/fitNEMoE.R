###################################################################
# function to fitting NEMoE with different initialization
###################################################################
#' @useDynLib NEMoE
#' @importFrom Rcpp sourceCpp
.fitNEMoE0 = function(X_list, Z, y, K, lambda1, alpha1 = 1, lambda2, alpha2 = 1,
                      gamma_init = NULL, beta_init = NULL, beta_max = 10,
                      EM_alg = "EM", init = "kmeans", itmax = 1e2, itmin = 3,
                      adapt = FALSE, btr = TRUE, stop_all = TRUE,
                      verbose = TRUE, early_stop = FALSE){


  n <- nrow(Z)
  L = length(X_list)
  q = ncol(Z)
  X_seg = .list2Matc(X_list)
  seg = X_seg$seg
  X = X_seg$X
  p_list = sapply(X_list,ncol)
  P = ncol(X)
  EM_list = 0:4
  names(EM_list) = c("EM","CEM", "SEM", "SAEM", "GEM")
  lambda1_mat = .paramExt(lambda1, seg, K)
  lambda1_list = .mat2Listr(lambda1_mat, seg)
  lambda2_mat = as.numeric(.paramExt(lambda2, q, 1))
  alpha1_mat = .paramExt(alpha1, L, K)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  if(is.null(names(X_list))){
    level_names = paste("L",1:L,sep = "")
    names(X_list) = level_names
  }else{
    level_names = names(X_list)
  }

  if(length(y) == n){
    y <- matrix(rep(y, L), ncol = L)
  }

  EM_opt = EM_list[EM_alg]

  #Initialize gamma
  if(is.null(gamma_init)){
    if(init == "rand"){
      a = rep(1, K)
      r_i = rDirichlet(n, a)
    }else{
      temp <- stats::kmeans(Z, centers = K)$cluster
      r_i <- .onehot(temp)
    }
  }else{
    r_i = calcProb(Z_1, gamma_init)
  }

  if(is.null(gamma_init)){
    gamma_init = sMulti(Z, r_i, lambda2_mat*L, rep(1,n), alpha2, beta_max)
  }

  #initialize beta
  if(is.null(beta_init)){
    beta_init = list()
    if(init == "glmnet"){
      for (i in 1:L) {
        beta_temp = sMulti(X_list[[i]], cbind(y[,i], 1-y[,i]),
                        rowMeans(lambda1_list[[i]]),
                        rep(1,n),mean(alpha1_mat[i,]), beta_max)
        beta_temp = 2*beta_temp[,2]
        beta_temp = matrix(rep(beta_temp, K), ncol = K)
        beta_init[[i]] = beta_temp
      }
    }else{
      for(i in 1:L){
        beta_temp <- matrix(0, nrow = (p_list[i] + 1), ncol = K)
        for(j in 1:K){
          beta_temp[,j] = 2*sMulti(X_list[[i]], cbind(y[,i],1-y[,i]),
                                lambda1_list[[i]][,j],
                                r_i[,j],alpha1_mat[i,j], beta_max)[,2]
        }
        beta_init[[i]] = beta_temp
      }
    }
  }
  beta_init0 = .list2Matr(beta_init)
  beta_seg = beta_init0$seg
  beta_init0 = beta_init0$X

  NEMoE_result = suppressWarnings(fitNEMoE0(X, seg, Z, y, K, lambda1_mat * K,
                                            lambda2_mat*L, alpha1_mat, alpha2,
                                            gamma_init, beta_init0, beta_max, EM_opt,
                                            itmax, itmin, adapt, btr,
                                            stop_all, verbose, early_stop))
  beta_fit1 = NEMoE_result[1:nrow(beta_init0),]
  beta_fit = .mat2Listr(beta_fit1, beta_seg)
  names(beta_fit) = level_names
  for(i in 1:L){
    if(is.null(colnames(X_list[[i]]))){
      temp_name = paste("gamma",1:ncol(X_list[[i]]), sep = "")
      rownames(beta_fit[[i]]) = c("intercept", temp_name)
    }else{
      rownames(beta_fit[[i]]) = c("intercept",
                              colnames(X_list[[i]]))
    }
  }

  gamma_fit = NEMoE_result[(nrow(beta_init0) + 1): (nrow(beta_init0) + nrow(gamma_init)), ]
  if(is.null(colnames(Z))){
    temp_name = paste("gamma",1:ncol(Z), sep = "")
    rownames(gamma_fit) = c("intercept", temp_name)
  }else{
    rownames(gamma_fit) = c("intercept",
                        colnames(Z))
  }

  r_i = array(0,dim = c(n, K, L))
  for(i in 1:L){
    r_i[,,i] = NEMoE_result[(nrow(beta_init0) + nrow(gamma_init) + n*(i-1) + 1):
                              (nrow(beta_init0) + nrow(gamma_init) + n*i),]
  }

  pen_fac1 = NEMoE_result[(nrow(beta_init0) + nrow(gamma_init) + n*L + 1) :
                           nrow(NEMoE_result),]
  pen_fac = .mat2Listr(pen_fac1, seg)

  LL_temp = calcLL(X, Z, y, seg, beta_fit1, gamma_fit, lambda1_mat*pen_fac1* K,
                   lambda2_mat*L, alpha1_mat, alpha2, T)
  colnames(LL_temp) <- c("LL_obs", "PLL_obs", "LL_complete", "PLL_complete")
  rownames(LL_temp) <- c(level_names, "All")
  if(L == 1){
    LL_temp = LL_temp[1,,drop = F]
  }

  return(list(beta = beta_fit, gamma = gamma_fit, r_i = r_i, pen_fac = pen_fac,
              LL = LL_temp))
}

###################################################################
# function to fitting NEMoE with different initialization
###################################################################

.fitNEMoE <- function(NEMoE, beta_init = NULL, gamma_init = NULL,
                      ...){

  X_list = NEMoE@Microbiome
  Z = NEMoE@Nutrition
  y = NEMoE@Response

  NEMoE_result = do.call(.fitNEMoE0,
                         c(list(X_list = X_list, Z= Z, y = y, K = NEMoE@K,
                                beta_init = beta_init, gamma_init = gamma_init),
                           NEMoE@params))

  NEMoE@NEMoE_output = NEMoE_result

  return(NEMoE)
}

###################################################################
# function to fitting NEMoE used in restart
###################################################################

.restart_NEMoE <- function(NEMoE){

  NEMoE@params$verbose = F
  NEMoE = .fitNEMoE(NEMoE)
  restart_result = NEMoE@NEMoE_output
  LL_temp =  restart_result$LL[nrow(restart_result$LL),2]
  return(list(beta = restart_result$beta,
              gamma = restart_result$gamma,
              PLL = LL_temp))

}

###################################################################
# function to fitting NEMoE with small EM
###################################################################
#' Fit NEMoE model
#' @description This function fit NEMoE model using EM algorithm.
#' @param NEMoE a NEMoE object contain data and parameters for fitting.
#' @param gamma_init initial value of parameters in gating network.
#' @param beta_init initial values of parameters in experts network.
#' @param BPPARAM A \code{BiocParallelParam} class object from.
#' @param num_restart Number of restart times.
#' @param restart_it Maximial number of iterations during restart.
#' the \code{BiocParallel} package is used. Default is SerialParam().
#' @param ... other parameters can be passed to fitNEMoE.
#' See \code{\link{createParameterList}}
#' @examples
#' data(NEMoE_example)
#' NEMoE_example
#' NEMoE_example <- fitNEMoE(NEMoE_example, num_restart = 0)
#' @return a NEMoE object contain fitted result including parameters
#' beta a list contain coefficients of experts network in all levels.
#' gamma a matrix of coefficients in gating network.
#' r_i an array of inherit partition result in EM algorithm.
#' pen_fac estimated penalty factor for lambda1 (adapt = T). The (balanced)
#' lambda1 in NEMoE is lambda1*pen_fac
#' latent estimated latent classes.
#'
#' @seealso \code{\link{NEMoE_buildFromList}},
#' \code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{createParameterList}}
#' @export

fitNEMoE = function(NEMoE, beta_init = NULL, gamma_init = NULL,
                    num_restart = 5, restart_it = 10,
                    BPPARAM = BiocParallel::SerialParam(),
                    ...){


  K <- NEMoE@K
  if(!(is.null(beta_init) & is.null(gamma_init))){

    return(.fitNEMoE(NEMoE, beta_init = beta_init, gamma_init = gamma_init))

  }else if(num_restart){

    km_flag = FALSE

    km_res0 <- stats::kmeans(NEMoE@Nutrition, centers = K)$cluster
    km_ari <- c()
    for(i in 1:num_restart){
      km_res <- stats::kmeans(NEMoE@Nutrition, centers = K)$cluster
      km_ari[i] <- mclust::adjustedRandIndex(km_res0, km_res)
    }

    if(min(km_ari) < 0.9){
      km_flag = TRUE
    }

    NEMoE_restart = NEMoE
    if(km_flag){
      NEMoE_restart@params$init = "kmeans"
    }else{
      NEMoE_restart@params$init = "rand"
    }
    NEMoE_restart@params$itmax = restart_it
    NEMoE_list = list()
    for(i in 1:num_restart){
      NEMoE_list[[i]] = NEMoE_restart
    }

    restart_res <- BiocParallel::bplapply(NEMoE_list, .restart_NEMoE,
                                          BPPARAM = BPPARAM)

    NEMoE_restart@params$init = "kmeans"
    restart_res[[num_restart + 1]] = .restart_NEMoE(NEMoE_restart)

    NEMoE_restart@params$init = "glmnet"
    restart_res[[num_restart + 2]] = .restart_NEMoE(NEMoE_restart)

    PLL_res = sapply(restart_res, `[[`, 3)

    idx = which.max(PLL_res)
    if(idx == (num_restart + 2)){
      NEMoE@params$init = "glmnet"
    }else if(idx == (num_restart + 1)){
      NEMoE@params$init = "kmeans"
    }else{
      NEMoE@params$init = "rand"
    }

    beta_init = restart_res[[idx]]$beta
    gamma_init = restart_res[[idx]]$gamma

    NEMoE <- .fitNEMoE(NEMoE, beta_init = beta_init, gamma_init = gamma_init)

  }else{

    NEMoE <- .fitNEMoE(NEMoE, beta_init = beta_init, gamma_init = gamma_init)

  }

  if(NEMoE@standardize){
    NEMoE <- .scale_back(NEMoE)
  }

  # Traceback unfiltered coefficient
  NEMoE <- .trace_filt(NEMoE)

  return(NEMoE)
}
