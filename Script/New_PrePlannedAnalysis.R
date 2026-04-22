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
d$choice_c   <- factor(d$choice, levels = c("skip", "help", "self"))
d$choice_ord <- ordered(d$choice, levels = c("skip", "help", "self"))

########################################################################
# RQ1: Math behaviours are ordered or categorical?
########################################################################
fit_cat <- brm(
  choice_c ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id),
  data = d, family = categorical(link = "logit"), cores = 4, warmup = 500, iter = 1000
  #backend = "cmdstanr", threads = threading(2)
)
fit_ord <- brm(
  choice_ord ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id),
  data = d, family = cumulative(link = "logit"), cores = 4, warmup = 500, iter = 1000
  #backend = "cmdstanr", threads = threading(2)
)

loo_cat <- loo(fit_cat); loo_ord <- loo(fit_ord)
comparison_rq1 <- loo_compare(loo_cat, loo_ord)
waic_cat <- waic(fit_cat); waic_ord <- waic(fit_ord)
comparison_rq1_waic <- loo_compare(waic_cat, waic_ord)

model_type <- ifelse(rownames(comparison_rq1)[1] == "fit_ord", "ordinal", "categorical")

########################################################################
# RQ2: State and trait measures as predictor of math behaviours
########################################################################

# Preliminary steps: checking pre-condition for 
# the inclusion of buoyancy and avoidance in the model
# 1. Checking fit indices factorial model for buoyancy and avoidance
# (single items were not simulated in example dataset,
# so items for both scales are simulated here for demonstration of analytical methods)
x = rnorm(1e3,0,1)

buoy1 = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2.5, -1.5, 0, 1, 1.5, 3, Inf), labels = c(1,2,3,4,5,6,7))
buoy2 = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -3, -2, -0.5, 0.5, 2, 2.5, Inf), labels = c(1,2,3,4,5,6,7))
buoy3 = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -1.5, 0.5, 0, 1.5, 2.5, 3, Inf), labels = c(1,2,3,4,5,6,7))
buoy4 = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2, -1, 0.5, 1.5, 2, 3.5, Inf), labels = c(1,2,3,4,5,6,7))

d_buoy = data.frame(buoy1, buoy2, buoy3, buoy4)


av1  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -1.5, 0, 1, Inf), labels = c(1,2,3,4))
av2  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2.5, 0.5, 1.5, Inf), labels = c(1,2,3,4))
av3  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -3, -0.5, 2, Inf), labels = c(1,2,3,4))
av4  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -1, 0.5, 1, Inf), labels = c(1,2,3,4))
av5  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -0.5, 0.5, 2.5, Inf), labels = c(1,2,3,4))
av6  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2, -1.5, 0, Inf), labels = c(1,2,3,4))
av7  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2, -1, 1.5, Inf), labels = c(1,2,3,4))
av8  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -1, 0, 1.5, Inf), labels = c(1,2,3,4))
av9  = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -1.5, 1, 2.5, Inf), labels = c(1,2,3,4))
av10 = cut(x+rnorm(1e3,0,1), breaks = c(-Inf, -2.5, -0.5, 1.5, Inf), labels = c(1,2,3,4))

d_av = data.frame(av1, av2, av3, av4, av5, av6, av7, av8, av9, av10)




# fitting models
model_buoy = 'buoyancy  =~ buoy1 + buoy2 + buoy3 + buoy4'
fitBuoy = cfa(model_buoy, data = d_buoy, ordered = names(d_buoy))
#summary(fitBuoy, fit.measures = T, standardized = T)

model_av = 'avoidance =~ av1  + av2  + av3  + av4 + av5 + av6 + av7 + av8 + av9 + av10'
fitAv = cfa(model_av, data = d_av, ordered = names(d_av))

# looking at fit indices including decisional rules
(fm = fitMeasures(fitBuoy, fit.measures=c("CFI","RMSEA")))
if(fm["cfi"]>0.95 & fm["rmsea"]<0.08){print("Fit indices ok: proceed!"); FitBuoy_ok = T
}else{print("Stop there!")}

(fm = fitMeasures(fitAv, fit.measures=c("CFI","RMSEA")))
if(fm["cfi"]>0.95 & fm["rmsea"]<0.08){print("Fit indices ok: proceed!"); FitAv_ok = T
}else{print("Stop there!")}

# 2. Checking reliability of latent factor
(rel = semTools::reliability(fitBuoy))
if(rel["alpha.ord",] > 0.70){print("Reliability ok: proceed!"); RelBuoy_ok = T
}else {print("Stop there!")}

(rel = semTools::reliability(fitAv))
if(rel["alpha.ord",] > 0.70){print("Reliability ok: proceed!"); RelAv_ok = T
}else {print("Stop there!")}


# predittori diretti condizionali al measurement model
buoy_incl = if (FitBuoy_ok & RelBuoy_ok) "+ buoyancy" else ""
avoid_incl = if (FitAv_ok & RelAv_ok) "+ avoidance" else ""


# path "a" — trait -> state
bf_worry   <- bf(worry ~ math_anxiety + math_ability_T0 + Gender + (1|id))
bf_scstate <- bf(sc_state ~ sc_trait + math_ability_T0 + Gender + (1|id))


# full mediation — starting model
if (model_type == "ordinal") {
  formula_full <- as.formula(paste(
    "choice_ord ~ worry + sc_state", buoy_incl, avoid_incl,
    "+ math_ability_T0 + Gender + (1|id)"
  ))
  bf_choice_full <- bf(formula_full, family = cumulative(link = "logit"))
} else {
  formula_full <- as.formula(paste(
    "choice_c ~ worry + sc_state", buoy_incl, avoid_incl,
    "+ math_ability_T0 + Gender + (1|id)"
  ))
  bf_choice_full <- bf(formula_full, family = categorical(link = "logit"))
}



# partial mediation — adding of direct paths trait → choice
if (model_type == "ordinal") {
  formula_partial <- as.formula(paste(
    "choice_ord ~ worry + sc_state + math_anxiety + sc_trait", buoy_incl, avoid_incl,
    "+ math_ability_T0 + Gender + (1|id)"
  ))
  bf_choice_partial <- bf(formula_partial, family = cumulative(link = "logit"))
} else {
  formula_partial <- as.formula(paste(
    "choice_c ~ worry + sc_state + math_anxiety + sc_trait", buoy_incl, avoid_incl,
    "+ math_ability_T0 + Gender + (1|id)"
  ))
  bf_choice_partial <- bf(formula_partial, family = categorical(link = "logit"))
}


# creation of model partially mediated vs. model fully mediated
fit_partial <- brm(
  bf_worry + bf_scstate + bf_choice_partial + set_rescor(FALSE),
  data = d, family = list(gaussian(), gaussian(), NULL),
  cores = 4, warmup = 500, iter = 1000
  #backend = "cmdstanr",
  #threads = threading(2)
)

fit_full <- brm(
  bf_worry + bf_scstate + bf_choice_full + set_rescor(FALSE),
  data = d, family = list(gaussian(), gaussian(), NULL),
  cores = 4, warmup = 500, iter = 1000
  #backend = "cmdstanr",
  #threads = threading(2)
)

# model comparison
loo_full <- loo(fit_full)
loo_partial <- loo(fit_partial)
comparison_mediation <- loo_compare(loo_full, loo_partial)

# model selection
partial_is_best <- rownames(comparison_mediation)[1] == "loo_partial"
models_equivalent  <- abs(comparison_mediation[2, "elpd_diff"]) < 
  2 * comparison_mediation[2, "se_diff"]

if (partial_is_best & !models_equivalent) {
  best_choice_model <- fit_partial
  mediation_type <- "partial"
} else {
  best_choice_model <- fit_full
  mediation_type <- "full"
}

cat("Mediation type selected:", mediation_type, "\n")

########################################################################
# indirect effects (path a*b)
########################################################################
#draws_med <- as_draws_df(fit_mediators)
draws_choice <- as_draws_df(best_choice_model)

# Path "a"
a1 <- draws_choice$b_worry_math_anxiety
a2 <- draws_choice$b_scstate_sc_trait

# Path "b"
if (model_type == "ordinal") {
  b1 <- draws_choice$b_worry
  b2 <- draws_choice$b_sc_state
} else {
  b1 <- draws_choice$b_muself_worry
  b2 <- draws_choice$b_muself_sc_state
}

indirect_MA <- a1 * b1
indirect_SC <- a2 * b2



# indirect effect 
results_indirect <- data.frame(
  Mediation            = c("MA tratto → worry → choice",
                           "SC tratto → sc_state → choice"),
  Indirect_Effect_Mean = c(mean(indirect_MA),            mean(indirect_SC)),
  CI_lower             = c(quantile(indirect_MA, 0.025), quantile(indirect_SC, 0.025)),
  CI_upper             = c(quantile(indirect_MA, 0.975), quantile(indirect_SC, 0.975)),
  PostProb_direction   = c(mean(indirect_MA < 0),        mean(indirect_SC > 0)),
  Model_type           = model_type,
  Mediation_type       = mediation_type
)

print(results_indirect)


########################################################################
# testing of the final model
########################################################################

if (model_type == "ordinal") {
  h_worry    <- hypothesis(best_choice_model, "worry < 0")
  h_sc_state <- hypothesis(best_choice_model, "sc_state > 0")
  
  hp_list <- list(
    list(predictor = "worry", direction = "<0", h = h_worry),
    list(predictor = "sc_state", direction = ">0", h = h_sc_state)
  )
  
  if (FitBuoy_ok & RelBuoy_ok) {
    h_buoyancy <- hypothesis(best_choice_model, "buoyancy > 0")
    hp_list <- c(hp_list, list(list(predictor = "buoyancy", direction = ">0", h = h_buoyancy)))
  }
  
  if (FitAv_ok & RelAv_ok) {
    h_avoidance <- hypothesis(best_choice_model, "avoidance < 0")
    hp_list<- c(hp_list, list(list(predictor = "avoidance", direction = "<0", h = h_avoidance)))
  }
  
  if (mediation_type == "partial") {
    h_ma_direct  <- hypothesis(best_choice_model, "math_anxiety < 0")
    h_sct_direct <- hypothesis(best_choice_model, "sc_trait > 0")
    hp_list <- c(hp_list,
                  list(list(predictor = "math_anxiety", direction = "<0", h = h_ma_direct)),
                  list(list(predictor = "sc_trait", direction = ">0", h = h_sct_direct)))
  }
  
  results_rq2 <- data.frame(
    predictor = sapply(hp_list, function(x) x$predictor),
    direction = sapply(hp_list, function(x) x$direction),
    Estimate  = sapply(hp_list, function(x) x$h$hypothesis$Estimate),
    CI.Lower  = sapply(hp_list, function(x) x$h$hypothesis$CI.Lower),
    CI.Upper  = sapply(hp_list, function(x) x$h$hypothesis$CI.Upper),
    Post.Prob = sapply(hp_list, function(x) x$h$hypothesis$Post.Prob)
  )
  
  # fdr correction for the whole set of predictors
  results_rq2$p_one_tailed <- 1 - results_rq2$Post.Prob
  results_rq2$p_fdr <- p.adjust(results_rq2$p_one_tailed, method = "fdr")
  
  print(results_rq2)
  
} else {
  
  h_worry_self    <- hypothesis(best_choice_model, "muself_worry < 0")
  h_sc_state_self <- hypothesis(best_choice_model, "muself_sc_state > 0")
  h_worry_help    <- hypothesis(best_choice_model, "muhelp_worry < 0")
  h_sc_state_help <- hypothesis(best_choice_model, "muhelp_sc_state > 0")
  
  hp_self <- list(
    list(predictor = "worry", direction = "<0", h = h_worry_self),
    list(predictor = "sc_state", direction = ">0", h = h_sc_state_self)
  )
  hp_help <- list(
    list(predictor = "worry", direction = "<0", h = h_worry_help),
    list(predictor = "sc_state", direction = ">0", h = h_sc_state_help)
  )
  
  if (FitBuoy_ok & RelBuoy_ok) {
    h_buoyancy_self <- hypothesis(best_choice_model, "muself_buoyancy > 0")
    h_buoyancy_help <- hypothesis(best_choice_model, "muhelp_buoyancy > 0")
    hp_self <- c(hp_self, list(list(predictor = "buoyancy", direction = ">0", h = h_buoyancy_self)))
    hp_help <- c(hp_help, list(list(predictor = "buoyancy", direction = ">0", h = h_buoyancy_help)))
  }
  
  if (FitAv_ok & RelAv_ok) {
    h_avoidance_self <- hypothesis(best_choice_model, "muself_avoidance < 0")
    h_avoidance_help <- hypothesis(best_choice_model, "muhelp_avoidance < 0")
    hp_self <- c(hp_self, list(list(predictor = "avoidance", direction = "<0", h = h_avoidance_self)))
    hp_help <- c(hp_help, list(list(predictor = "avoidance", direction = "<0", h = h_avoidance_help)))
  }
  
  if (mediation_type == "partial") {
    h_ma_self  <- hypothesis(best_choice_model, "muself_math_anxiety < 0")
    h_sct_self <- hypothesis(best_choice_model, "muself_sc_trait > 0")
    h_ma_help  <- hypothesis(best_choice_model, "muhelp_math_anxiety < 0")
    h_sct_help <- hypothesis(best_choice_model, "muhelp_sc_trait > 0")
    hp_self <- c(hp_self,
                  list(list(predictor = "math_anxiety", direction = "<0", h = h_ma_self)),
                  list(list(predictor = "sc_trait",     direction = ">0", h = h_sct_self)))
    hp_help <- c(hp_help,
                  list(list(predictor = "math_anxiety", direction = "<0", h = h_ma_help)),
                  list(list(predictor = "sc_trait",     direction = ">0", h = h_sct_help)))
  }
  
  # function to creating dataset for contrasts
  make_df <- function(hp_list, contrast_label) {
    data.frame(
      predictor = sapply(hp_list, function(x) x$predictor),
      contrast  = contrast_label,
      direction = sapply(hp_list, function(x) x$direction),
      Estimate  = sapply(hp_list, function(x) x$h$hypothesis$Estimate),
      CI.Lower  = sapply(hp_list, function(x) x$h$hypothesis$CI.Lower),
      CI.Upper  = sapply(hp_list, function(x) x$h$hypothesis$CI.Upper),
      Post.Prob = sapply(hp_list, function(x) x$h$hypothesis$Post.Prob)
    )
  }
  
  results_rq2 <- rbind(make_df(hp_self, "self_vs_skip"),
                       make_df(hp_help, "help_vs_skip"))
  
  results_rq2$p_one_tailed <- 1 - results_rq2$Post.Prob
  
  # FDR correction applied separately within each contrast
  for (ctr in unique(results_rq2$contrast)) {
    idx <- results_rq2$contrast == ctr
    results_rq2$p_fdr[idx] <- p.adjust(results_rq2$p_one_tailed[idx], method = "fdr")
  }
  
  print(results_rq2)
  
  # secondary constrasts self vs help
  h_diff <- hypothesis(best_choice_model,
                       c("muself_worry < muhelp_worry",
                         "muself_sc_state > muhelp_sc_state"))
  print(h_diff)
}








