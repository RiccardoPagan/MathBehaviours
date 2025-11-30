# empty workspace
rm(list=ls()); gc()

#load required packages
library(brms)
library(dplyr)
library(lme4)
library(effects)
#library(ordinal)

Path <- getwd()

# import dataset
d <- read.csv("Data/SimulatedDatasetBehaviours.csv", sep = ",")


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
#my_formula_cat <- choice_c ~ math_anxiety + (1|id)
d$choice_c  <- factor(d$choice,
                      levels = c("skip", "help", "self"))

my_formula_ord <- choice_ord ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id)
#my_formula_ord <- choice_ord ~ math_anxiety + (1|id)
d$choice_ord <- ordered(d$choice, 
                        levels = c("skip", "help", "self"))

# model with behaviours as categorical variable
fit_cat <- brm(
  formula = my_formula_cat,
  data = d,
  family = categorical(link = "logit"),
  cores = 4
)

# model with choice as ordered (skip < help < self)
fit_ord <- brm(
  formula = my_formula_ord,
  data = d,
  family = cumulative(link = "logit"),
  cores = 4
)


# model comparison
loo_cat <- loo(fit_cat)
loo_ord <- loo(fit_ord)
comparison1 <- loo_compare(loo_cat, loo_ord)

waic_cat <- waic(fit_cat)
waic_ord <- waic(fit_ord)
comparison2 <- loo_compare(waic_cat, waic_ord)


# output results
setwd(paste0(Path, "/outputs"))

sink("output.txt")
cat("=========================================\n")
cat("Categorical model\n")
cat("=========================================\n\n")
summary(fit_cat)

cat("\n\n")

cat("=========================================\n")
cat("Ordered model\n")
cat("=========================================\n\n")
summary(fit_ord)
sink()
closeAllConnections()


sink("LOO.txt")
cat("=========================================\n")
cat("LOO categorical model\n")
cat("=========================================\n\n")
loo_cat

cat("\n\n")

cat("=========================================\n")
cat("LOO ordered model\n")
cat("=========================================\n\n")
loo_ord

sink()
closeAllConnections()


sink("WAIC.txt")
cat("=========================================\n")
cat("WAIC categorical model\n")
cat("=========================================\n\n")
waic_cat

cat("\n\n")

cat("=========================================\n")
cat("WAIC ordered model\n")
cat("=========================================\n\n")
waic_ord

sink()
closeAllConnections()



sink("comparison.txt")
cat("=========================================\n")
cat("LOO COMPARISON\n")
cat("=========================================\n\n")
comparison1

cat("\n\n")

cat("=========================================\n")
cat("WAIC COMPARISON\n")
cat("=========================================\n\n")
comparison2

sink()
closeAllConnections()

########################################################################
# RQ2: State and trait measures as predictor of math behaviours
########################################################################

# model with behaviours as categorical variable (apparently best model)
fit_cat <- brm(
  formula = my_formula_cat,
  data = d,
  family = categorical(link = "logit"),
  cores = 4
)

summary(fit_cat)
plot(fit_cat)







