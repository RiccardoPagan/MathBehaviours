# empty workspace
rm(list=ls())
gc()

#load required packages
library(brms)
library(dplyr)
library(lme4)
library(effects)
#library(ordinal)

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

my_formula <- choice ~ buoyancy + avoidance + worry + pc_state + math_anxiety + pc_trait + math_ability + Gender + (1|id)

# model with behaviours as categorical variable
fit_cat <- brm(
  formula = my_formula,
  data = d,
  family = categorical(link = "logit"),
  cores = 4
)


# model with choice as ordered (skip < help < self)
d$choice <- factor(d$choice, 
                   levels = c("skip", "help", "self"), 
                   ordered = TRUE)

# model with behaviours as ordered variable
fit_ord <- brm(
  formula = my_formula,
  data = d,
  family = cumulative(link = "logit"),
  cores = 4
)



# model comparison
comparison1 <- loo(fit_ord, fit_cat)

loo_cat = loo(fit_cat)
loo_ord = loo(fit_ord)

comparison2 <- loo_compare(loo_cat, loo_ord)















#model with choice as categorical variable

fit1 <- brm(choice ~ math_anxiety + buoyancy + avoidance + pc + worry +
              metacognition + math_ability + (1|ID), data=d,
            family=categorical(), cores=4)

summary(fit1)
plot(fit1)


########################################################################
# RQ2: Math behaviours as predictor of math performance
########################################################################

d <- subset(d, session == 2)
d_ns <- subset(d, choice != "skip")

# model without skip
fit2 <- glmer(
  accuracy ~ choice + math_ability + (1|ID),
  data = d_ns,
  family = binomial(link = "logit")
)

summary(fit2)
plot(allEffects(fit2),multiline=T,ci.style="auto")


# bayesan model to try to stabilize the skip coefficient with a prior (not the ideal solution)
fit_bay <- brm(
  accuracy ~ choice + math_ability + (1 | ID),
  data = d,
  family = bernoulli(link = "logit"),
  prior = prior(normal(0, 1), class = "b"), 
  cores = 4, chains = 4
)

summary(fit_bay)
plot(fit_bay)






