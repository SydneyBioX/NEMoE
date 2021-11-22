#' Transformation of phyloseq data
#' @param ps a phyloseq object.
#' @param tax_list Get taxonomic levels of data. By default, return the taxonomic of Phylum, Order, Family, Genus and ASV/OTU.
#' @param prop_trans a logic variable to determine whether the data are as proportion. By default is TRUE.
#' @return A list of data of different taxonomic levels.
#' @examples
#' data(PD)
#' ps1 <- psGather(PD$ps)
#' @importFrom stats na.omit
#' @export

psGather <- function(ps, tax_list = c("Phylum","Order","Family","Genus","ASV"),
                     prop_trans = T){

  tax_rank <- phyloseq::rank_names(ps)

  idx <- match(tax_list, tax_rank)

  n_tax <- na.omit(idx)

  data_list <- list()

  for(i in 1:length(n_tax)){

    ps_temp <- phyloseq::tax_glom(ps, taxrank = tax_rank[idx[i]])

    if(prop_trans){

      ps_temp <- phyloseq::transform_sample_counts(ps_temp,
                                                   function(x){x / sum(x)})

    }

    X_temp <- ps_temp@otu_table

    name_temp <- ps_temp@tax_table[,idx[i]]

    colnames(X_temp) <- name_temp

    data_list[[i]] <- X_temp

    names(data_list)[i] <- tax_rank[idx[i]]
  }

  if(any(c("ASV","OTU") %in% tax_list)){

    if(prop_trans){
      ps_temp <- phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})
    }else{
      ps_temp <- ps
    }

    gen_name <- ps_temp@tax_table[,6]

    name_ASV <- paste( gen_name,
                       paste("ASV.", 1:length(gen_name), sep = ""),
                       sep = ":")

    X_temp = ps_temp@otu_table
    colnames(X_temp) <- name_ASV

    data_list[["ASV"]] = X_temp
  }

  return(data_list)

}

#' Transformation of compositional data
#' @param X a matrix to transform.
#' @param eps when elements of X smaller than eps, will be eps to avoid invalid value in log transformation.
#' @param method The method of transformation.
#' If method = "asin" will perform arcsine transformation.
#' If method = "clr" will perform central log transformation.
#' If method = "comp" will use composition data.
#' @param scale a logical variable to determine whether scale the transformed data.
#' @return A matrix of transformed data.
#' @export

compTransform <- function(X, eps = 1e-4, method = "comp", scale = F){

  if(method == "clr"){

    X_log <- log(X + 0.5)

    X_trans <- X_log - rowMeans(X_log)

  }else if(method == "asin"){

    X_comp <- X/rowSums(X)

    X_trans <- asin(sqrt(X_comp))
  }else if(method == "none"){

    X_trans <- X

  }else{

    X_trans <- X/rowSums(X)
  }

  if(scale){
    X_trans <- scale(X_trans)
  }

  return(X_trans)
}

#' Filter the composition matrix
#' @param X the matrix input to filter
#' @param thresh_func a function to calculate statistics for filtering.
#'  By default is variance.
#' @param thresh If the statistics of variable smaller than threshold,
#' the corresponding variable will filter out. By default is 1e-4.
#' @param id_out Whether output the id of variable kept in final output.
#' @return A list of filtered data and/or index.
#' @export
#' @importFrom stats var
filterComp <- function(X, thresh_func = var, thresh = 1e-4, id_out = FALSE){

  X <- as.matrix(X)

  X_stat <- apply(X, 2, thresh_func)

  id <- X_stat > thresh
  names(id) <- colnames(X)

  X <- X[, (X_stat > thresh)]

  if(id_out){
    return(list(X = X, id = id))
  }else{
    return(list(X = X))
  }
}

# Calculate how many variables included in the matrix
# If version = "all", return how many non-zero rows
# Else return how many non-zero elements

.calcdf <- function(X, version = "all", eps = 1e-8, intercept = TRUE){

  X <- as.matrix(X)

  if(intercept){
    p <- (nrow(X) - 1)
    X <- X[2:(p + 1),,drop = F]
  }else{
    p <- nrow(X)
  }

  if(version == "all"){

    var_tf <- (abs(X) > eps)

    var_tf <- apply(var_tf, 1, any)

    return(sum(as.numeric(var_tf)))

  }else{

    var_tf <- (abs(X) > eps)
    return(sum(apply(var_tf, 2, as.numeric)))
  }

}

#' Calculate degree of freedom
#' @description This function calculate degree of freedom of NEMoE object.
#' @param NEMoE A NEMoE object with fitted result.
#' @param output whether degree of freedom with combined
#'  components or compute within each component.
#' @return a list of degree of freedom in gating network and experts network.
#' @examples
#' data(NEMoE_example)
#' calcdf(NEMoE_example)
#' @export

calcdf <- function(NEMoE, output = "all"){

  K <- NEMoE@K
  L <- length(NEMoE@Microbiome)
  gamma <- NEMoE@NEMoE_output$gamma
  beta <- NEMoE@NEMoE_output$beta

  df_gamma <- .calcdf(gamma)

  if(output == "all"){
    df_beta <- matrix(0, nrow = L, ncol = 1)
    for(i in 1:L){
      df_beta[i,] <- .calcdf(beta[[i]])
    }
  }else{
    df_beta <- matrix(0, nrow = L, ncol = K)
    for(i in 1:L){
      for(j in 1:K){
        df_beta[i,j] <- .calcdf(beta[[i]][,j,drop = F], "comp")
      }
    }
  }
  return(list(df_gamma = df_gamma, df_beta = df_beta))


}

# softmax function
.softmax <- function(x){

  return(exp(x)/sum(exp(x)))

}

# Clipping function to make the value within a range
.clipping <- function(x, max_x, min_x){

  x[x > max_x] = max_x

  x[x < min_x] = min_x

  return(x)
}

# Logistic function
.logistic <- function(x){

  return(1 / (1 + exp(-x)))

}

# Derivitive of logistic function
.d_logistic <- function(x, beta){

  return(x - x * .logistic(x %*% beta))

}

# A soft threshold function (Solution for lasso)
.soft_thresh <- function(x, lambda){

  y <- x

  y[x>lambda] <- x - lambda

  y[x<(-lambda)] <- x + lambda

  y[(x<= lambda)&(x >= -lambda)] = 0

  return(y)
}

# One hot encoding to convert a factor result to one hot matrix

.onehot <- function(x){

  K <- max(x)
  n <- length(x)

  one_hot <- matrix(0, nrow = n, ncol = K)

  for(i in 1:n){
    one_hot[i,x[i]] <- 1
  }

  return(one_hot)
}

#' Convert label to factor of matrix
#' @description This function convert a soft probability matrix to factor or one hot matrix.
#' @param prob a matrix to convert.
#' @param idx the type of output. If idx= TRUE will
#'  return the index of largest probability of each row. If idx= FALSE will return a one-hot matrix.
#' @return a index of largest probability or one hot matrix.
#' @export

cvtLabel <- function(prob, idx = FALSE){

  prob<- as.matrix(prob)
  n <- nrow(prob)
  K <- ncol(prob)

  prob_arg <- apply(prob, 1, which.max)

  if(idx){
    return(prob_arg)
  }else{

    prob_thresh <- matrix(0, nrow = n, ncol = K)

    for(i in 1:n){
      prob_thresh[i,prob_arg[i]] <- 1
    }

    return(prob_thresh)
  }
}

# A function for converting a list to matrix put into function NEMoEfit0.

.list2Matc <- function(X_list){

  L <- length(X_list)
  n <- nrow(X_list[[1]])
  seg = c()
  X <- matrix(0, nrow = n, ncol = 1)
  for(i in 1:L){
    seg[i] = ncol(X_list[[i]])
    X <- cbind(X, X_list[[i]])
  }
  X <- X[,-1]
  return(list(X = X, seg= seg))
}

.list2Matr <- function(X_list){

  L <- length(X_list)
  p <- ncol(X_list[[1]])
  seg = c()
  X <- matrix(0, nrow = 1, ncol = p)
  for(i in 1:L){
    seg[i] = nrow(X_list[[i]])
    X <- rbind(X, X_list[[i]])
  }
  X <- X[-1,]
  return(list(X = X, seg= seg))
}

.mat2Listc <- function(X, seg){
  L <- length(seg)

  X_list <- list()
  for(i in 1:L){
    X_list[[i]] <- X[,1:seg[i]]
    X <- X[,-c(1:seg[i])]
  }
  return(X_list)
}

.mat2Listr <- function(X, seg){
  L <- length(seg)

  X_list <- list()
  for(i in 1:L){
    X_list[[i]] <- X[1:seg[i],]
    X <- X[-c(1:seg[i]),]
  }
  return(X_list)
}

.paramExt <- function(param, p_list, K){

  P <- sum(p_list)
  L <- length(p_list)

  if(length(param) == 1){
    param_mat = matrix(param, nrow = P, ncol = K)
  }else if(length(param) == K){
    param_mat = matrix(rep(param, P), nrow = P, ncol = K, byrow = T)
  }else if(length(param) == L){
    param_mat <- matrix(0, nrow = 1, ncol = K)
    for(i in 1:L){
      param_temp <- matrix(rep(param[i], p_list[i]*K), ncol = K)
      param_mat <- rbind(param_mat, param_temp)
    }
    param_mat = param_mat[-1,]
  }else if(length(param) == L*K){
    param_mat <- matrix(0, nrow = 1, ncol = K)
    for(i in 1:L){
      param_temp <- matrix(rep(param[i, ], p_list[i]), ncol = K, byrow = T)
      param_mat <- rbind(param_mat, param_temp)
    }
    param_mat = param_mat[-1,]
  }else{
    temp <- ceiling((P*K)/length(param))
    param_temp <- rep(param, temp)[1:(P*K)]
    param_mat = matrix(param_temp, ncol = K)
  }
  return(param_mat)
}

#' @importFrom stats sd
.colSds <- function(X){
  return(apply(X,2,sd))
}

.mathcNEMoE <- function(Z, Z_new){

  Z_match <- matrix(0, nrow = nrow(Z_new), ncol = ncol(Z))

  if(is.null(colnames(Z))){
    Z_match[,1:ncol(Z_new)] = Z_new
  }else if(is.null(colnames(Z_new))){
    Z_match[,1:ncol(Z_new)] = Z_new
  }else{
    Z_int <- intersect(colnames(Z_new), colnames(Z))
    idx_match <- match(Z_int, colnames(Z))

    Z_match[,na.omit(idx_match)] <- Z_new
  }

  return(Z_match)
}

.sample_data <- function(ps){
  return(phyloseq::sample_data(ps))
}

.scale_back <- function(NEMoE){

  L <- length(NEMoE@Microbiome)

  .transformation <- NEMoE@.transformation
  NEMoE_result <- NEMoE@NEMoE_output
  beta = NEMoE_result$beta
  gamma = NEMoE_result$gamma
  Z_sd <- .transformation$sd_Z
  Z_mu <- .transformation$mu_Z
  gamma[2:nrow(gamma),] = gamma[2:nrow(gamma),]/Z_sd
  gamma[1,] = gamma[1,] - colSums((Z_mu * gamma[2:nrow(gamma), ]))

  X_mu <- .transformation$mu_X
  X_sd <- .transformation$sd_X

  for(i in 1:L){
    beta_temp = beta[[i]]
    beta_temp[2:nrow(beta_temp),] = beta_temp[2:nrow(beta_temp),]/X_sd[[i]]
    beta_temp[1,] = beta_temp[1,] -
      colSums((X_mu[[i]]*beta_temp[2:nrow(beta_temp),]))
    beta[[i]] = beta_temp
  }
  NEMoE_result$beta = beta
  NEMoE_result$gamma = gamma
  NEMoE@NEMoE_output <- NEMoE_result

  return(NEMoE)
}

.trace_filt <- function(NEMoE){

  L <- length(NEMoE@Microbiome)
  K <- NEMoE@K

  beta <- list()

  id0 <- NEMoE@.transformation$keepid
  p_list <- sapply(id0, length)
  beta0 <- NEMoE@NEMoE_output$beta
 for(i in 1:L){

    idx_temp <- unname(c(TRUE,id0[[i]]))
    names_temp <- names(id0[[i]])
    if(is.null(names_temp)){
      names_temp <- paste("V",idx_temp, sep = "")
    }
    names_temp <- c("intercept",names_temp)

    beta_temp <- matrix(0, nrow = (p_list[i] + 1), ncol = K)
    beta_temp[as.logical(idx_temp),] = beta0[[i]]
    rownames(beta_temp) <- names_temp
    beta[[i]] <- beta_temp
 }
  names(beta) <- names(beta0)
  NEMoE@NEMoE_output$beta = beta
  return(NEMoE)
}
