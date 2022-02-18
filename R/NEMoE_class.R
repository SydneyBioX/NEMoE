#' check validity
#' @param object A NEMoE object.
check_validity <- function(object){

  n_M = sapply(object@Microbiome, nrow)
  n_N = nrow(object@Nutrition)
  n_R = length(object@Response)

  if(!all(n_M == n_M[1])){
    errors <- c("Numbers of samples in microbiome data are not equal!")
    return(errors)
  }
  if(!all(n_M == n_N)){
    errors <- c("Numbers of samples in microbiome data not equal
                to nutirtion data!")
    return(errors)
  }
  if(!all(n_M == n_R)){
    errors <- c("Numbers of samples in microbiome data not equal
                to heath outcome!")
    return(errors)
  }
  if(!is.null(object@taxLevel)){
    if(length(object@Microbiome) != length(object@taxLevel)){
      errors <- c("Selected taxa level should be the same as
                  names in microbiome list!")
      return(errors)
    }
  }
  if(is.integer(object@K)){
    errors = c("Number of latent classes should be integer!")
  }

  return(TRUE)
}

#' @importClassesFrom S4Vectors DataFrame DFrame
setClassUnion("data.frameORNULL", c("data.frame", "NULL"))

setClassUnion("characterORNULL", c("character", "NULL"))

setClassUnion("listORNULL", c("list", "NULL"))

#' NEMoE class
#' @slot Microbiome A list of microbiome matrix, rows are samples.
#' @slot Nutrition A matrix of nutrition intake.
#' @slot Response A vector of health outcome.
#' @slot params A list of parameters for fitting NEMoE.
#' @slot K A number of fitting componenets in NEMoE (>=2).
#' @slot cvParams A list of parameters for selecting parameters of NEMoE
#'  using cross validation.
#' @slot cvResult A dataframe of cross validation results of different parameters
#'  of fitting NEMoE.
#' @slot NEMoE_output A list NEMoE fitting result.
#' @slot taxLevel A character indicates selected taxa level.
#' @slot ResponseLevel A character indicates levels of response variables.
#' @slot standardize A Logical variable indicate whether standardize input data.
#' @importFrom methods new
NEMoE_class <- setClass("NEMoE",
                   slots = c(Microbiome = "list",
                             Nutrition = "matrix",
                             Response = "vector",
                             params = "listORNULL",
                             K = "numeric",
                             cvParams = "listORNULL",
                             cvResult = "listORNULL",
                             NEMoE_output = "list",
                             taxLevel = "characterORNULL",
                             taxTab = "data.frameORNULL",
                             ResponseLevel = "characterORNULL",
                             standardize = "logical",
                             .transformation = "listORNULL"
                             ),
                   prototype = list(
                     params = list(),
                     cvParams = list(),
                     .transformation = list(method = "none")
                   ))

#' @importFrom S4Vectors coolcat
#'

setMethod("show", "NEMoE", function(object) {

  L = length(object@Microbiome)
  Microbiome_char = "Microbiome number of variables: "
  if(is.null(object@taxLevel)){
    for(i in 1:L){
      Microbiome_char = paste0(Microbiome_char,
                               dim(object@Microbiome[[i]])[2],
                               ", ")
    }
    Microbiome_char = paste0(Microbiome_char, "\n")
  }else{
    for(i in 1:L){
      Microbiome_char = paste0(Microbiome_char,
                               dim(object@Microbiome[[i]])[2],
                               "(", object@taxLevel[i],")",
                               ", ")
    }
    Microbiome_char = paste0(Microbiome_char, "\n")
  }

  cat("Numer of samples: ", nrow(object@Microbiome[[1]]), '\n')
  cat(Microbiome_char)
  cat("Nutrition dim:", dim(object@Nutrition), "\n")
  cat("Number of latent class K: ", object@K, '\n')
  S4Vectors::coolcat("NEMoE output names (%d): %s\n",
                     names(object@NEMoE_output))

  S4Vectors::coolcat("param names (%d): %s\n", names(object@params))
  S4Vectors::coolcat("cross validation param names (%d): %s\n",
                     names(object@cvParams))
})

#' @importFrom S4Vectors setValidity2
#'
S4Vectors::setValidity2("NEMoE", check_validity)

