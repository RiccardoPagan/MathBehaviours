rm(list=ls()); gc()

library(brms)
library(bayestestR)
library(lavaan)
library(dplyr)
library(semTools)

Path <- getwd()
d <- read.csv("SimulatedDatasetBehaviours.csv", sep = ",")

d$choice <- paste0(d$skip, d$help, d$self)
d$choice <- case_when(
  d$choice == "100" ~ "skip",
  d$choice == "010" ~ "help",
  d$choice == "001" ~ "self",
  TRUE ~ NA
)
d <- subset(d, session == 2)
#d$choice_c   <- factor(d$choice, levels = c("skip", "help", "self"))
d$choice_ord <- ordered(d$choice, levels = c("skip", "help", "self"))


library(brms)
library(tidyverse)

fit_ord_base <- brm(
  choice_ord ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1 | id),
  data    = d,
  family  = cumulative(link = "logit"),
  iter = 1000, warmup = 500,
  cores = 4
)


fit_ord_1 <- brm(
  choice_ord ~ buoyancy + avoidance + cs(worry) + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1 | id),
  data    = d,
  family  = cumulative(link = "logit"),
  iter = 1000, warmup = 500,
  cores = 4
)


fit_ord_base <- add_criterion(fit_ord_base, c("loo", "waic"))
fit_ord_1 <- add_criterion(fit_ord_1, c("loo", "waic"))

loo_compare <- loo_compare(
  fit_ord_base,
  fit_ord_1
)
print(loo_compare, simplify = FALSE)
