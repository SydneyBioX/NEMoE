#' @title 16S gut microbiome Parkinson's disease dataset
#' @description  A list of two type of data (A phyloseq object and a list that can be put into NEMoE.
#' Microbiome data：
#' The phyloseq object contain 995 taxa and 168 samples.
#' The data list contain 168 samples with 4 levels(5, 16, 24, 36, 101) variables respectively.
#' Nutrition data: 168 samples and 27 variables.
#' PD state: PD = 1 and HC = 0.
#' @format A list of two type of data:
#' \describe{
#'   \item{ps}{A phyloseq object with 168 samples and 995 taxa}
#'   \item{data_list}{A list of Microbiome, Nutrition and response variables}
#' }
#' @usage data(PD, package = 'NEMoE')
"PD"

#' @title An example of fitted NEMoE object
#' @description An example of fitted NEMoE object
#' Generate from following code
#' data(PD)
#' Microbiome = PD$data_list$Microbiome
#' Nutrition = PD$data_list$Nutrition
#' Response = PD$data_list$Response
#' NEMoE = NEMoE_buildFromList(Microbiome, Nutrition, Response,
#' lambda1 = c(0.005, 0.014, 0.016, 0.023, 0.025),lambda2 = 0.02,
#' alpha1 = 0.5, alpha2 = 1, adapt = T, btr = T, stop_all = T, itmax = 1e3,
#' verbose = T, beta_max = 1e2, early_stop = F)
#' NEMoE = fitNEMoE(NEMoE, num_restart = 0)
#' @usage data(NEMoE_example, package = 'NEMoE')
"NEMoE_example"
