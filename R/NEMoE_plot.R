###################################################################
# Plot fitting result of gating network
###################################################################
#' Plot gating network
#' @description This function plot the PCA of fitted latent class
#' and their corresponding loadings.
#' @param NEMoE_obj a NEMoE object with fitted output.
#' @param PCs Visualization of selected Principal components.
#' @return A graph of PCA plot of nutrition intake and estimated latent classes
#' and its corresponding loadings.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic coord_flip
#' @examples
#' data(NEMoE_example)
#' plotGating(NEMoE_example)
#' @export

plotGating <- function(NEMoE_obj, PCs = c(1,2)){
  if(!length(NEMoE_obj@NEMoE_output)){
    error <- "Please fitting NEMoE before plot."
    return(TRUE)
  }

  gamma <- NEMoE_obj@NEMoE_output$gamma
  y <- as.factor(NEMoE_obj@Response)

  colnames(gamma) <- paste("latent",1:ncol(gamma), sep = "")
  gamma_df <- as.data.frame(gamma)
  gamma_df <- gamma_df[-1,]
  gamma_df$var_name = rownames(gamma_df)
  gamma_df_melt = reshape2::melt(gamma_df, id.vars = "var_name")
  colnames(gamma_df_melt)[2] = "latent"

  Z <- NEMoE_obj@Nutrition
  n <- nrow(Z)
  Z1 <- cbind(rep(1,n), Z)
  latent <- calcProb(Z1, gamma)
  latent <- as.factor(cvtLabel(latent, idx = TRUE))

  Z_pca <- stats::prcomp(Z, scale. = T)
  Z_pca <- as.data.frame(Z_pca$x[,PCs])
  Z_pca$latent <- latent
  p1 <- ggplot(Z_pca) +
    geom_point(aes(x = Z_pca[,1], y = Z_pca[,2],
                   color = .data$latent, shape = y))+
    theme_classic() + labs(x = paste0("PC",PCs[1]), y = paste0("PC",PCs[2]))
  p2 <- ggplot(gamma_df_melt) +
    geom_bar(aes(x = .data$var_name, y = .data$value, fill = .data$latent),
             stat = "identity") + coord_flip() +
    labs(x= "", y ="") + theme_classic()

  return(list(p1, p2))

}
###################################################################
# Plot fitting result of experts network
###################################################################
#' Plot experts network
#' @description This function plot the estimated coefficients of each latent class
#' on each level.
#' @param NEMoE_obj a NEMoE object with fitted output.
#' @return A list of graph of fitted coefficients of each latent class
#' on each level.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic
#' @examples
#' data(NEMoE_example)
#' plotExperts(NEMoE_example)

#' @export

plotExperts <- function(NEMoE_obj){

  beta = NEMoE_obj@NEMoE_output$beta
  L <- length(beta)
  p_list = list()
  for(i in 1:L){
    beta_temp <- beta[[i]]
    beta_temp_df <- as.data.frame(beta_temp)
    colnames(beta_temp_df) <- paste("latent",1:ncol(beta_temp_df), sep = "")
    beta_temp_df <- beta_temp_df[-1,]
    beta_temp_df$var_name = rownames(beta_temp_df)
    beta_temp_df_melt = reshape2::melt(beta_temp_df, id.vars = "var_name")
    colnames(beta_temp_df_melt)[2] = "latent"

    p_temp <-  ggplot(beta_temp_df_melt) +
      geom_bar(aes(x = .data$var_name, y = .data$value, fill = .data$latent),
               stat = "identity", position = "dodge") +
      labs(x= "", y ="") + theme_classic()
    p_list[[i]] = p_temp
  }

  return(p_list)
}
