# empty workspace and RAM memory
rm(list=ls()); gc()

#load required packages
if (!"brms" %in% installed.packages()) install.packages("brms")
if (!"bayestestR" %in% installed.packages()) install.packages("bayestestR")
library(brms)
library(bayestestR)
library(dplyr)
library(lme4)
library(effects)
library(ordinal)

Path <- getwd()

# import dataset
d <- read.csv("SimulatedDatasetBehaviours.csv", sep = ",")


# data preparation
d$choice = paste0(d$skip,d$help,d$self)
d$choice = case_when(
  d$choice=="100" ~ "skip",
  d$choice=="010" ~ "help",
  d$choice=="001" ~ "self",
  TRUE ~ NA
)

########################################################################
# RQ1: Math behaviours are ordered or categorical?
########################################################################
d <- subset(d, session == 2)

my_formula_cat <- choice_c ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id)
d$choice_c  <- factor(d$choice,
                      levels = c("skip", "help", "self"))

my_formula_ord <- choice_ord ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id)
d$choice_ord <- ordered(d$choice, 
                        levels = c("skip", "help", "self"))

# model with behaviours as categorical variable
fit_cat <- brm(
  formula = my_formula_cat,
  data = d,
  family = categorical(link = "logit"),
  cores = 7
)

# model with choice as ordered (skip < help < self)
fit_ord <- brm(
  formula = my_formula_ord,
  data = d,
  family = cumulative(link = "logit"),
  cores = 7
)


# model comparison
loo_cat <- loo(fit_cat)
loo_ord <- loo(fit_ord)
comparison1 <- loo_compare(loo_cat, loo_ord)

waic_cat <- waic(fit_cat)
waic_ord <- waic(fit_ord)
comparison2 <- loo_compare(waic_cat, waic_ord)


# results output in a .txt file
# setwd(paste0(Path, "/outputs"))
# 
# sink("output.txt")
# cat("=========================================\n")
# cat("Categorical model\n")
# cat("=========================================\n\n")
# summary(fit_cat)
# 
# cat("\n\n")
# 
# cat("=========================================\n")
# cat("Ordered model\n")
# cat("=========================================\n\n")
# summary(fit_ord)
# sink()
# closeAllConnections()
# 
# 
# sink("LOO.txt")
# cat("=========================================\n")
# cat("LOO categorical model\n")
# cat("=========================================\n\n")
# loo_cat
# 
# cat("\n\n")
# 
# cat("=========================================\n")
# cat("LOO ordered model\n")
# cat("=========================================\n\n")
# loo_ord
# 
# sink()
# closeAllConnections()
# 
# 
# sink("WAIC.txt")
# cat("=========================================\n")
# cat("WAIC categorical model\n")
# cat("=========================================\n\n")
# waic_cat
# 
# cat("\n\n")
# 
# cat("=========================================\n")
# cat("WAIC ordered model\n")
# cat("=========================================\n\n")
# waic_ord
# 
# sink()
# closeAllConnections()
# 
# 
# 
# sink("comparison.txt")
# cat("=========================================\n")
# cat("LOO COMPARISON\n")
# cat("=========================================\n\n")
# comparison1
# 
# cat("\n\n")
# 
# cat("=========================================\n")
# cat("WAIC COMPARISON\n")
# cat("=========================================\n\n")
# comparison2
# 
# sink()
# closeAllConnections()

########################################################################
# RQ2: State and trait measures as predictor of math behaviours
########################################################################

# select the best model based on model comparison
best_model <- fit_cat  # or fit_cat based on comparison results
model_type <- "categorical" # or 'categorical' based on comparison results 

if (model_type == "ordinal") {
  
  # for ordinal model there is a single coefficient per predictor
  h1 <- hypothesis(best_model, "buoyancy > 0")
  h2 <- hypothesis(best_model, "avoidance < 0")
  h3 <- hypothesis(best_model, "worry < 0")
  h4 <- hypothesis(best_model, "sc_state > 0")
  
} else if (model_type == "categorical") {
  
  # for categorical model there are separate coefficients for each comparison
  # test for "self" vs "skip" comparison (self-reliance)
  h1_self <- hypothesis(best_model, "muself_buoyancy > 0")
  h2_self <- hypothesis(best_model, "muself_avoidance < 0")
  h3_self <- hypothesis(best_model, "muself_worry < 0")
  h4_self <- hypothesis(best_model, "muself_sc_state > 0")
  
  # Test for "help" vs "skip" comparison
  h1_help <- hypothesis(best_model, "muhelp_buoyancy > 0")
  h2_help <- hypothesis(best_model, "muhelp_avoidance < 0")
  h3_help <- hypothesis(best_model, "muhelp_worry < 0")
  h4_help <- hypothesis(best_model, "muhelp_sc_state > 0")
  
  # test differences between self and help (exploratory contrasts)
  h_diff <- hypothesis(best_model, 
                       c("muself_buoyancy > muhelp_buoyancy",
                         "muself_sc_state > muhelp_sc_state",
                         "muself_avoidance < muhelp_avoidance"))
  

  results_rq2 <- data.frame(
    predictor  = rep(c("buoyancy", "avoidance", "worry", "sc_state"), each = 2),
    contrast   = rep(c("self_vs_skip", "help_vs_skip"), times = 4),
    direction  = rep(c(">0", "<0", "<0", ">0"), each = 2),
    Estimate   = c(h1_self$hypothesis$Estimate,  h1_help$hypothesis$Estimate,
                   h2_self$hypothesis$Estimate,  h2_help$hypothesis$Estimate,
                   h3_self$hypothesis$Estimate,  h3_help$hypothesis$Estimate,
                   h4_self$hypothesis$Estimate,  h4_help$hypothesis$Estimate),
    CI.Lower   = c(h1_self$hypothesis$CI.Lower,  h1_help$hypothesis$CI.Lower,
                   h2_self$hypothesis$CI.Lower,  h2_help$hypothesis$CI.Lower,
                   h3_self$hypothesis$CI.Lower,  h3_help$hypothesis$CI.Lower,
                   h4_self$hypothesis$CI.Lower,  h4_help$hypothesis$CI.Lower),
    CI.Upper   = c(h1_self$hypothesis$CI.Upper,  h1_help$hypothesis$CI.Upper,
                   h2_self$hypothesis$CI.Upper,  h2_help$hypothesis$CI.Upper,
                   h3_self$hypothesis$CI.Upper,  h3_help$hypothesis$CI.Upper,
                   h4_self$hypothesis$CI.Upper,  h4_help$hypothesis$CI.Upper),
    Post.Prob  = c(h1_self$hypothesis$Post.Prob, h1_help$hypothesis$Post.Prob,
                   h2_self$hypothesis$Post.Prob, h2_help$hypothesis$Post.Prob,
                   h3_self$hypothesis$Post.Prob, h3_help$hypothesis$Post.Prob,
                   h4_self$hypothesis$Post.Prob, h4_help$hypothesis$Post.Prob)
  )
  
  results_rq2$p_one_tailed <- 1 - results_rq2$Post.Prob
  
  # FDR correction applied separately within each contrast set
  for (ctr in unique(results_rq2$contrast)) {
    idx <- results_rq2$contrast == ctr
    results_rq2$p_fdr[idx] <- p.adjust(results_rq2$p_one_tailed[idx], method = "fdr")
  }
  
  print(results_rq2)
  cat("\n--- Exploratory: Self vs. Help contrasts ---\n")
  print(h_diff)
}

