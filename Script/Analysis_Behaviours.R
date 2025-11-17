# empty workspace
rm(list=ls())
gc()

#load required packages
library(brms)
library(dplyr)
library(lme4)

# import dataset
d <- read.csv("Data/SimulatedDatasetBehaviours.csv", sep = ",")

# data preparation
library(dplyr)
d$choice = paste0(d$skip,d$help,d$self)
d$choice = case_when(
  d$choice=="100" ~ "skip",
  d$choice=="010" ~ "help",
  d$choice=="001" ~ "self",
  TRUE ~ NA
)

########################################################################
# RQ1: State and trait measures as predectors of math behaviours
########################################################################

fit1 <- brm(choice ~ math_anxiety + buoyancy + avoidance + pc + worry +
            metacognition + math_ability + (1|ID), data=d,
          family=categorical(), cores=4)

########################################################################
# RQ2: State and trait measures as predectors of math behaviours
########################################################################

d <- subset(d, session == 2)

fit2 <- glmer(
  accuracy ~ choice + math_ability + (1|ID),
  data = d,
  family = binomial(link = "logit")
)



