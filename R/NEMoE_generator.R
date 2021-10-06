###################################################################
# function to generate simulation data
###################################################################
#' Generate mixture distribution data
#' @description This function generate a data follows a mixture distribution with different
#' sample size and variables in the dataset.
#'
#' @param n Number of samples when generating the dataset. By default is 200.
#' @param p Number of variables for experts network input. By default is 30.
#' @param q Number of variables for gating netwrok By defult is 20.
#' @param s1 Number of non-zeros coefficient in
#'  experts network input. By default is 5.
#' @param s2 Number of non-zeros coefficient in
#'  gating network input By default is 5.
#' @param p_L A numeric vector of length (L-1), each entries indicate
#'  number of variables in each level.
#' @param K0 Number of components for latent class in dataset. By default is 2.
#' @param Sigma Covariance matrix for gating network input.
#'  If it is NULL, will take identity matrix as covariance.
#' @param eta Coefficient of separation parameters in
#'  generating data for gating networks. By default is 0.5.
#' @param c_g Coefficient of signal strength parameters in
#'  generating data for gating networks. By default is 1.
#' @param c_e Coefficient of signal strength parameters in
#'  experts data for gating networks. By default is 1.
#' @param link the method for generating response \code{y}.
#' If link = "probit", use mixture of probit model.
#' If version  = "logit", use mixture of logistic model. By default is logit.
#' @param method The transformation method used for
#'  construct relationship in experts network.
#' If method = "comp", use prepositional data.
#'  If method = "asin", use arcsin transformed compositional data.
#' If method = "clr", use central log ratio transformed compositional data.
#'  By default is "comp".
#' @param scale Logical variable to indicate whether to use
#'  scaled coefficient. By default is TRUE
#' @param gen_Micro A character indicates which model used in generate
#' microbiome data, can be chosen from "zinLDA", "dm" and "mgauss", means
#' zero-inflated latent Dirichelet allocation model, Dirichlet multinomial model
#' and multivariate gaussian model.
#' @param prev_filt The threshold of prevalence of selected that
#'  have non-zero coefficients. By default is 0.3.
#' @param var_filt The threshold of variance of selected that have
#'  non-zero coefficients. By default is 1e-6.
#' @param beta_max Maximal number of coefficients for experts network.
#' @param fix_X Fixed microbiome input matrix. If NULL,
#'  will generate using zinLDA model.
#' @param ... other parameters can be passed to genNEMoE.
#' i.e. parameters in zinLDA (K = 5, Alpha = 10, Pi = 0.4, a = 0.05, b = 10)
#' @return A list contain the generated microbiome data\code{X},
#'  nutrition data\code{W},
#' health response \code{y}, coefficients of experts network \code{beta},
#'  coefficients of gating network \code{gamma},
#' simulated observed logits \code{pi}, simulated latent group \code{latent}
#'  and simulated response probability \code{y_prob}.
#' @examples
#' dat <- genNEMoE(n = 10, p = 10000, q = 30)
#' @importFrom stats runif var
#' @export
#'
genNEMoE <- function(n = NULL, p = NULL, q = 30, K0 = 2, Sigma = NULL,
                     eta = 0.5, c_g = 1, c_e = 1, s1 = 3, s2 = 4, p_L = c(10,20,50),
                     fix_X = NULL, gen_Micro = "zinLDA", prev_filt = 0.3, var_filt = 1e-6,
                     method = "comp", scale = T, link = "probit", beta_max = 100, ...){


  c_g = c_g*2
  if(length(fix_X)){
    X = fix_X
  }else{
    N <- runif(n, 5000, 25000)
    X = .genMicro(N = N, n = n, p = p, gen_Micro, ...)
  }

  if(!length(Sigma)){
    Sigma = diag(rep(1, q))
  }

  n <- nrow(X)
  p <- ncol(X)

  Nutri_gen <- .genGating(n = n, q = q, K = K0, s_gamma = s1, Sigma = Sigma,
                          eta = eta, c_g = c_g)

  Z = Nutri_gen$Z
  beta_V <- Nutri_gen$beta_V
  tax_tab <- .genLevel(p_in = p, p_out = p_L)


  X_bin <- (X > 0)
  if(gen_Micro != "mgauss"){
    X_comp <- compTransform(X, method = "comp", scale = F)
  }else{
    X_comp <- X
  }

  prev_X <- colSums(X_bin)/n
  var_X <- apply(X_comp,2,var)

  filt_idx <- which((prev_X > prev_filt)&(var_X > var_filt))
  X_eff <- X[,filt_idx]
  tax_tab_eff <- tax_tab[filt_idx,]

  Experts_gen <- .genExperts(X = X_eff, tax_tab = tax_tab_eff,
                             K = K0, s = s2, c_e = c_e,
                             method = method, scale = scale,
                             beta_max = beta_max)

  X_list = Experts_gen$X_list
  beta_W = Experts_gen$beta_W

  pi <- .genLatent(Z, K0, beta_V)
  latent <- cvtLabel(pi)

  L = length(X_list)
  y <- matrix(0, nrow = n, ncol = L)
  y_prob <- matrix(0, nrow = n, ncol = L)

  for(l in 1:L){

    gen_y_list <- .genY(X_list[[l]], Z, pi = pi, beta_W = beta_W[[l]],
                        K = K0, link = link)

    latent <- gen_y_list$latent

    y[,l] <- gen_y_list$y

    y_prob[,l] <- gen_y_list$y_prob
  }

  y_obs <- round(rowMeans(y_prob))


  return(list(X_count = X, X_list = X_list, tax_tab = tax_tab,
              W = Z, gamma = beta_V, beta = beta_W, pi = pi,
              latent = latent, y = y_obs, y_prob = y_prob))

}


###################################################################
# function to generate nutrition data
###################################################################
#' @importFrom stats rmultinom
.genNutri <- function(n = 200, q = 30, K = 2,
                      pi = NULL, mu = NULL, Sigma = NULL){

  if(!length(pi)){
    pi = rep(1/K, K)
  }

  z_latent <- rmultinom(n, 1, pi)

  if(!length(Sigma)){
    Sigma = diag(rep(1, q))
  }

  if(!length(mu)){

    mu = MASS::mvrnorm(n = K, mu = rep(0,q), Sigma = diag(rep(1, q)))

  }

  X <- matrix(0, nrow = n, ncol = q)

  if(K > 1){
    for(i in 1:K){
      n_temp <- sum(z_latent[i,])
      X_temp <- MASS::mvrnorm(n = n_temp, mu = mu[i,], Sigma = Sigma)
      X[z_latent[i,] == 1,] = X_temp
    }
  }else if (K == 1){
    X = MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  }

  return(X)

}

###################################################################
# function to generate microbiome data
###################################################################
#' @useDynLib NEMoE
.genMicro_zinLDA <- function(N, n, p, K = 5, Alpha = 10, Pi = 0.4,
                      a = 0.05, b = 10){

  taxa <- paste0("Taxa ",seq(p))
  communities <- paste0("Community ", seq(K))
  samples <- paste0("Sample ", seq(n))

  Delta = do.call(cbind, lapply(1:p, function(i)
    {stats::rbinom(n = K, size = 1, prob = Pi)}))
  rownames(Delta) <- communities
  colnames(Delta) <- taxa

  deltaCheck = T
  while(deltaCheck){
    V.rs = which(colSums(Delta) == K)
    n.V.rs = length(V.rs)
    Delta[,V.rs] = do.call(cbind, lapply(1:n.V.rs, function(i)
      {stats::rbinom(n = K, size = 1, prob = Pi)}))
    deltaCheck = length(which(colSums(Delta) == K))
  }

  Beta <- ifelse(Delta == 1, 0, NA)

  for (j in 1:K) {
    Q_j_index <- which(Delta[j,] == 0)
    Q_j <- stats::rbeta(length(Q_j_index)-1, a, b)

    for (l in seq_along(Q_j_index[1:length(Q_j_index)-1])) {
      Beta[j, Q_j_index[l]] <- Q_j[l]*prod(1 - Q_j[0:(l-1)])
    }

    Beta[j, Q_j_index[length(Q_j_index)]] <- 1 - sum(Beta[j,], na.rm = T)
  }

  Theta <- rDirichlet(n = n, a = rep(Alpha, K))
  rownames(Theta) <- samples
  colnames(Theta) <- communities

  generateSample <- function(N, theta, beta){

    sample <- vector(length = N)
    z_d <- vector(length = N)

    for (n in 1:N) {
      # For each N in the sample
      z_n <- stats::rmultinom(1, 1, theta)
      w_n <- stats::rmultinom(1, 1, beta[which(z_n == 1),])

      z_d[n] <- which(z_n == 1)
      sample[n] <- colnames(beta)[which(w_n == 1)]
      names(z_d)[n] <- sample[n]
    }

    return(list("sample" = sample,
                "z" = z_d))
  }

  cohort <- vector(mode = "list", length = n)
  z <- vector(mode = "list", length = n)

    for (d in 1:n) {
    sample.d <- generateSample(N[d], Theta[d,], Beta)
    cohort[[d]] <- sample.d[["sample"]]
    z[[d]] <- sample.d[["z"]]
  }

  sampleTaxaFreq <- lapply(cohort, table)
  sampleTaxaMatrix <- matrix(data = 0, nrow = n, ncol = p)
  rownames(sampleTaxaMatrix) <- samples
  colnames(sampleTaxaMatrix) <- taxa

  for (d in 1:n) {
    sampleTaxaMatrix[d, names(sampleTaxaFreq[[d]])] <- sampleTaxaFreq[[d]]
  }

  return(sampleTaxaMatrix)

}

.genMicro_dm <- function(N, n, p, a){

  if (missing(N)){
    N <- runif(n, 5000, 25000)
  }

  if (length(N) == 1)
      N <- stats::rpois(n, N)

  if (missing(a)){
    a <- stats::rnorm(p, mean = 14, sd = 4)
  }else if(length(a) == 1){
    a <- rep(a, p)
  }else if(length(a) <p){
    a[(length(a)+1):p] = mean(a)
  }else{
    a <- a[1:p]
  }
  a <- abs(a)

  P <- rDirichlet(n, a * (1 - 0.3)/0.3)

  X <- matrix(0, n, p)

  for (i in 1:n) X[i, ] <- rmultinom(1, N[i], P[i, ])

  return(X)


}

.genMicro_mgauss <- function(n, p, a = 0, b = 0){

  mu <- rep(a, p)
  Sigma <- (1 - b) * diag(rep(1,p)) + b*matrix(1, nrow = p, ncol = p)

  X = MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  return(X)
}

.genMicro <- function(N, n, p, gen_Micro, ...){

  if(gen_Micro == "zinLDA"){
    X <- .genMicro_zinLDA(N, n, p, ...)
  }else if(gen_Micro == "dm"){
    X <- .genMicro_dm(N, n, p, ...)
  }else{
    X <- .genMicro_mgauss(n, p, ...)
  }

  return(X)
}

###################################################################
# function to generate taxa levels
###################################################################

.genLevel <- function(p_in = 80, p_out = c(10,30,50)){

  L <- length(p_out)
  L <- L + 1

  tax_tab <- matrix(0, nrow = p_in, ncol = L)
  if(L == 1){
    tax_tab[,1] = 1:p_in
  }else{
    tax_tab[,L] = 1:p_in
    for(l in seq(L-1, 1, -1)){

      tax_tab[,l] = sample(p_out[l],p_in, replace = T)

    }
  }
  return(tax_tab)
}

###################################################################
# function to generate nutrition class
###################################################################

.genLatent <- function(Z, K, beta_V, r){

  n <- nrow(Z)
  p1 <- ncol(Z)

  latent_odd <- Z %*% beta_V

  latent_prob <- t(apply(latent_odd, MARGIN = 1, .softmax))

  return(latent_prob)

}

###################################################################
# function to generate response in each submodel
###################################################################
#' @importFrom stats rnorm rbinom
.genSub <- function(X, beta_W, link = "prob"){

  n <- nrow(X)
  p <- ncol(X)

  logit_odd <- X %*% beta_W

  logit_prob <- exp(logit_odd) / (1 + exp(logit_odd))

  if(link == "logit"){

    binom <- purrr::partial(rbinom, n = 1, size = 1)
    y <- 1 - apply(logit_prob, MARGIN = 1, binom)

  }else{

    y <- round(logit_prob)

  }

  return(list(y = y, prob = logit_prob))

}

###################################################################
# function to generate response
###################################################################
#' @useDynLib NEMoE
.genY <- function(X, Z, beta_V = NULL, beta_W, K, pi = NULL, link = "probit"){

  n <- nrow(X)

  if(is.null(pi)){
    pi <- .genLatent(Z, K, beta_V)
  }


  if(link == "logit"){

    latent <- rSample(pi)

  }else{

    latent <- cvtLabel(pi)
  }

  y <- rep(0, n)
  y_prob <- rep(0, n)

  for(k in 1:K){

    nk <- sum(latent[,k] == 1)

    y_sub <- .genSub(X[latent[,k] == 1,], beta_W[,k])

    y_prob[latent[,k] == 1] <- y_sub$prob

    y[latent[,k] == 1] = y_sub$y

  }

  return(list(pi = pi, latent = latent, y = y, y_prob = y_prob))

}


###################################################################
# function to generator sparse index in each model
# using a "most" exclusive principle
###################################################################

.genIdx <- function(p_in, p_out, s){

  p <- length(p_in) + length(p_out)

  if(s > 1){
    s_real <- s
  }else{
    s_real <- ceiling(p*s)
  }

  if(!length(p_out)){
    idx <- sample(p_in, s_real)
  }else if (length(p_out) < s_real){
    idx <- sample(p_in, s_real - length(p_out))
    idx <- union(idx, p_out)
  }else{
    idx <- sample(p_out, s_real)
  }

  p_in <- union(p_in, idx)
  p_out <- setdiff(p_out, idx)

  return(list(p_in = p_in, p_out = p_out, idx = idx))
}

###################################################################
# function to generator mean of nutrition distribution
###################################################################
#' @importFrom stats rnorm
.genMu <- function(K, p1, idx, niu){

  if (K == 1){
    mu = rep(0, p1)
  }else{
    mu = rep(0, p1)
    temp = rnorm(idx)
    mu[idx] = sign(temp) * niu
  }

  return(mu)
}


###################################################################
# function to generator coefficients
###################################################################
#' @importFrom stats rnorm
.genBeta <- function(p, idx, c = 1){

  beta <- rep(0, p)

  beta_temp <- rnorm(idx)

  beta[idx] <- sign(beta_temp) * c

  return(beta)
}


###################################################################
# function to generate covariate and coefficients of gating network
###################################################################
.genGating <- function(n, q, K, s_gamma, eta, c_g, Sigma = NULL){

  if(!length(Sigma)){
    Sigma = diag(rep(1, q))
  }

  beta_V <- matrix(0, nrow = q, ncol = K)
  for(i in 1:K){
    beta_V[sample(1:q, s_gamma),i] <- c_g*sign(rnorm(s_gamma))
  }

  latent <- sample(1:K, n, replace = T)

  W <- matrix(0, nrow = n, ncol = q)

  for(i in 1:n){

    W[i,] = MASS::mvrnorm(1, mu = eta*beta_V[,latent[i]], Sigma)
  }

  mu_z <- eta*t(beta_V)

  Z <- .genNutri(n, q, K = K, mu = mu_z, Sigma = Sigma)

  return(list(Z = Z, beta_V = beta_V))
}

##################################################################################
# function to generate multi-level covariates and
# coefficients for experts network
##################################################################################
#' @importFrom stats sd
.genExperts <- function(X, tax_tab, K, s, c_e,
                        method = "comp", scale = F, beta_max = 20){

  n <- nrow(X)
  p <- ncol(X)

  tax_tab <- as.matrix(tax_tab)
  L <- ncol(tax_tab)

  beta_W <- list()
  X_list <- list()

  X_temp <- X
  colnames(X_temp) <- paste("L",L,"_",1:p, sep = "")
  X_temp <- compTransform(X_temp, method = method, scale = scale)
  X_list[[L]] <- X_temp
  sd_X <- apply(X_temp,2, sd)

  beta_base <- matrix(0, nrow = p, ncol = K)

  X_bin <- (X_temp > 0)
  prev_X <- colSums(X_bin)/nrow(X)
  var_X <- apply(X_temp, 2, var)
  sd_X <- apply(X_temp, 2, sd)
  idx0 <- which(prev_X > 0.3 & var_X > 1e-5)
  p_in1 <- setdiff(1:p, idx0)
  p_out1 <- idx0

  for(k in 1:K){

    idx_gen = .genIdx(p_in1, p_out1, s)

    beta_temp = .genBeta(p, idx_gen$idx, c_e)

    beta_base[,k] <- .clipping(1/sd_X * beta_temp, beta_max, -beta_max)

  }
  beta_W[[L]] <- beta_base

  if(L > 1){
    for(l in 1:(L - 1)){

      l_i <- sort(unique(tax_tab[,l]),decreasing = F)
      X_temp <- matrix(0, nrow = n, ncol = length(l_i))
      beta_temp <- matrix(0, nrow = length(l_i), ncol = K)

      for( i in 1:length(l_i)){
        X_temp[,i] <- rowSums(X[,which(tax_tab[,l] == l_i[i]), drop = F])
        beta_temp[i,] <- colMeans(beta_base[which(tax_tab[,l] == l_i[i])
                                            ,, drop = F])
      }

      colnames(X_temp) <- paste("L",l,"_",l_i, sep = "")
      X_temp <- compTransform(X_temp, method = method, scale = scale)
      X_list[[l]] <- X_temp
      beta_W[[l]] <- beta_temp
    }
  }else{
    X_list[[1]] <- X_temp
  }

  return(list(X_list = X_list, beta_W = beta_W))

}
