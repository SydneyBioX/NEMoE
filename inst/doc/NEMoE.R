## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ggplot2)
library(NEMoE)
library(phyloseq)

## -----------------------------------------------------------------------------
data("PD")

## -----------------------------------------------------------------------------
Microbiome = PD$data_list$Microbiome
Nutrition = PD$data_list$Nutrition
Response = PD$data_list$Response

## -----------------------------------------------------------------------------
NEMoE = NEMoE_buildFromList(Microbiome, Nutrition, Response, K = 2, 
                            lambda1 = c(0.005, 0.014, 0.016, 0.023, 0.025),
                            lambda2 = 0.02, alpha1 = 0.5, alpha2 = 0.5)

## -----------------------------------------------------------------------------
NEMoE = fitNEMoE(NEMoE)

## -----------------------------------------------------------------------------
getLL(NEMoE)

## -----------------------------------------------------------------------------
coef.Gating <- getCoef(NEMoE)$coef.gating

## -----------------------------------------------------------------------------
plotGating(NEMoE)

## -----------------------------------------------------------------------------
coef.Experts <- getCoef(NEMoE)$coef.experts

## -----------------------------------------------------------------------------
plotExperts(NEMoE)

## -----------------------------------------------------------------------------
Microbiome_gen <- Microbiome$Genus

## -----------------------------------------------------------------------------
NEMoE_gen = NEMoE_buildFromList(Microbiome_gen, Nutrition, Response, K = 2, 
                            lambda1 = 0.028, lambda2 = 0.02,
                            alpha1 = 0.5, alpha2 = 0.5)

## -----------------------------------------------------------------------------
NEMoE_gen <- fitNEMoE(NEMoE_gen, restart_it = 20)

## -----------------------------------------------------------------------------
plotGating(NEMoE_gen)

## -----------------------------------------------------------------------------
p <- plotExperts(NEMoE_gen)[[1]]
p + theme(axis.text.x = element_text(angle = 45))

## -----------------------------------------------------------------------------
sample_PD <- sample_data(PD$ps)
Response <- sample_PD$PD
Nutrition <- sample_PD[,-ncol(sample_PD)]

## -----------------------------------------------------------------------------
NEMoE_ps <- NEMoE_buildFromPhyloseq(ps = PD$ps, Nutrition = scale(Nutrition),
                                 Response = Response,  K = 2, 
                                 lambda1 = c(0.005, 0.014, 0.016, 0.023, 0.025),
                                 lambda2 = 0.02, alpha1 = 0.5, alpha2 = 0.5,
                                 filtParam = list(prev = 0.7, var = 1e-5),
                                 transParam = list(method = "asin", scale = T),
                                 taxLevel = c("Phylum","Order","Family",
                                              "Genus","ASV"))

## -----------------------------------------------------------------------------
NEMoE_ps <- fitNEMoE(NEMoE_ps)

## -----------------------------------------------------------------------------
plotGating(NEMoE_ps)

## -----------------------------------------------------------------------------
plotExperts(NEMoE_ps)

## -----------------------------------------------------------------------------
calcCriterion(NEMoE, "all")

## -----------------------------------------------------------------------------
NEMoE <- cvNEMoE(NEMoE)

## -----------------------------------------------------------------------------
lambda1_choose <- NEMoE@cvResult$lambda1_choose
lambda2_choose <- NEMoE@cvResult$lambda2_choose

## -----------------------------------------------------------------------------
NEMoE <- setParam(NEMoE, lambda1 = lambda1_choose, lambda2 = lambda2_choose)

## -----------------------------------------------------------------------------
NEMoE <- fitNEMoE(NEMoE)

## -----------------------------------------------------------------------------
plotGating(NEMoE)

## -----------------------------------------------------------------------------
plotExperts(NEMoE)

## -----------------------------------------------------------------------------
sessionInfo()

