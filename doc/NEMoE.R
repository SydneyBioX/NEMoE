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
                            lambda2 = 0.02, alpha1 = 0.5, alpha2 = 0.5,
                            cvParams = createCVList(g1 = 10, shrink = 0.4,
                                                    track = F))

## -----------------------------------------------------------------------------
NEMoE = fitNEMoE(NEMoE)

## -----------------------------------------------------------------------------
getLL(NEMoE)

## -----------------------------------------------------------------------------
coef.Gating <- getCoef(NEMoE)$coef.gating

## -----------------------------------------------------------------------------
p_list <- plotGating(NEMoE)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]])

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]])

## -----------------------------------------------------------------------------
coef.Experts <- getCoef(NEMoE)$coef.experts

## -----------------------------------------------------------------------------
p_list <- plotExperts(NEMoE)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]] + ggtitle("Coefficients of Phylum level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]] + ggtitle("Coefficients of Order level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[3]] + ggtitle("Coefficients of Family level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[4]] + ggtitle("Coefficients of Genus level") + theme(axis.text.x = element_text(angle = 45, size = 5))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[5]] + ggtitle("Coefficients of ASV level") + theme(axis.text.x = element_blank())) 

## -----------------------------------------------------------------------------
Microbiome_gen <- Microbiome$Genus

## -----------------------------------------------------------------------------
NEMoE_gen = NEMoE_buildFromList(Microbiome_gen, Nutrition, Response, K = 2, 
                            lambda1 = 0.028, lambda2 = 0.02,
                            alpha1 = 0.5, alpha2 = 0.5)

## -----------------------------------------------------------------------------
NEMoE_gen <- fitNEMoE(NEMoE_gen, restart_it = 20)

## -----------------------------------------------------------------------------
p_list <- plotGating(NEMoE_gen)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]])

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]])

## ----fig.height=8, fig.width=8------------------------------------------------
p <- plotExperts(NEMoE_gen)[[1]]
p + theme(axis.text.x = element_text(angle = 45, size = 5))

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
p_list <- plotGating(NEMoE_ps)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]])

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]])

## ----fig.height=8, fig.width=8------------------------------------------------
p_list <- plotExperts(NEMoE)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]] + ggtitle("Coefficients of Phylum level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]] + ggtitle("Coefficients of Order level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[3]] + ggtitle("Coefficients of Family level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[4]] + ggtitle("Coefficients of Genus level") + theme(axis.text.x = element_text(angle = 45, size = 5))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[5]] + ggtitle("Coefficients of ASV level") + theme(axis.text.x = element_blank())) 

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
p_list <- plotGating(NEMoE_ps)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]])

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]])

## -----------------------------------------------------------------------------
p_list <- plotExperts(NEMoE)

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[1]] + ggtitle("Coefficients of Phylum level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[2]] + ggtitle("Coefficients of Order level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[3]] + ggtitle("Coefficients of Family level") + theme(axis.text.x = element_text(angle = 45))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[4]] + ggtitle("Coefficients of Genus level") + theme(axis.text.x = element_text(angle = 45, size = 5))) 

## ----fig.height=8, fig.width=8------------------------------------------------
print(p_list[[5]] + ggtitle("Coefficients of ASV level") + theme(axis.text.x = element_blank())) 

## -----------------------------------------------------------------------------
sessionInfo()

