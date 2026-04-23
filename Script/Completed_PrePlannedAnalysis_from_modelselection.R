rm(list = ls()); gc()

library(brms)
library(bayestestR)
library(lavaan)
library(dplyr)
library(semTools)
library(posterior)
library(tidyverse)

Path <- getwd()
d <- read.csv("SimulatedDatasetBehaviours.csv", sep = ",")

d$choice <- paste0(d$skip, d$help, d$self)
d$choice <- case_when(
  d$choice == "100" ~ "skip",
  d$choice == "010" ~ "help",
  d$choice == "001" ~ "self",
  TRUE ~ NA_character_
)
d <- subset(d, session == 2)
d$choice_ord <- ordered(d$choice, levels = c("skip", "help", "self"))

# optional standardization of continuous predictors
vars_to_scale <- intersect(
  c("buoyancy", "avoidance", "worry", "sc_state",
    "math_anxiety", "sc_trait", "math_ability_T0"),
  names(d)
)
d[vars_to_scale] <- lapply(d[vars_to_scale], function(x) as.numeric(scale(x)))

#########################################################################
# FIT HELPER
########################################################################

fit_choice_model <- function(formula, data) {
  brm(
    formula = formula,
    data    = data,
    family  = cumulative(link = "logit"),
    iter    = 1000,
    warmup  = 500,
    chains  = 4,
    cores   = 4,
    seed    = 220426
  )
}

########################################################################
# HELPERS
########################################################################

compare_by_elpd <- function(loo_obj) {
  best_name <- rownames(loo_obj)[1]
  ratio <- abs(loo_obj[1, "elpd_diff"] / loo_obj[1, "se_diff"])
  list(best = best_name, clear = is.finite(ratio) && ratio > 2, ratio = ratio)
}

summarise_draws <- function(x, direction = c("positive", "negative"), label = NA_character_) {
  direction <- match.arg(direction)
  data.frame(
    parameter = label,
    Estimate = mean(x),
    CI.Lower = unname(quantile(x, 0.025)),
    CI.Upper = unname(quantile(x, 0.975)),
    Post.Prob = if (direction == "positive") mean(x > 0) else mean(x < 0),
    stringsAsFactors = FALSE
  )
}

extract_univariate_choice_draws <- function(fit, predictor) {
  dr <- as_draws_df(fit)
  nm <- names(dr)

  pat_standard <- paste0("^b_", predictor, "$")
  pat_cs_1 <- paste0("^bcs_", predictor, "($|\\[)")
  pat_cs_2 <- paste0("^bcs_.*_", predictor, "($|\\[)")

  hits <- nm[grepl(pat_standard, nm) | grepl(pat_cs_1, nm) | grepl(pat_cs_2, nm)]

  if (length(hits) == 0) {
    stop(paste("No coefficient found for predictor:", predictor))
  }

  out <- lapply(hits, function(h) dr[[h]])
  names(out) <- hits
  out
}

fit_screening_model <- function(data_person, items, latent_name, alpha_cut = 0.50) {
  if (!all(items %in% names(data_person))) {
    return(list(
      fit_ok = FALSE,
      rel_ok = FALSE,
      included = FALSE,
      fit_measures = NA,
      alpha_ord = NA,
      note = paste("Item columns not found for", latent_name)
    ))
  }

  model_string <- paste0(latent_name, " =~ ", paste(items, collapse = " + "))
  fit <- cfa(model_string, data = data_person[, items, drop = FALSE], ordered = items)
  fm <- fitMeasures(fit, fit.measures = c("cfi", "rmsea"))

  rel <- tryCatch(
    semTools::reliability(fit),
    error = function(e) NULL
  )

  alpha_ord <- if (!is.null(rel) && "alpha.ord" %in% rownames(rel)) {
    as.numeric(rel["alpha.ord", 1])
  } else {
    NA_real_
  }

  fit_ok <- isTRUE(fm["cfi"] > 0.95 && fm["rmsea"] < 0.08)
  rel_ok <- isTRUE(!is.na(alpha_ord) && alpha_ord > alpha_cut)

  list(
    fit = fit,
    fit_ok = fit_ok,
    rel_ok = rel_ok,
    included = fit_ok && rel_ok,
    fit_measures = fm,
    alpha_ord = alpha_ord,
    note = NULL
  )
}

########################################################################
# PRELIMINARY SCREENING FOR BUOYANCY AND AVOIDANCE
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




# looking at fit indices including decisional rules (Initial fit check)
(fm = fitMeasures(fitBuoy, fit.measures=c("CFI","RMSEA")))
if(fm["cfi"] > 0.95 & fm["rmsea"] < 0.08){print("Fit indices ok: proceed!"); FitBuoy_ok_initial = TRUE
} else {print("Initial fit insufficient - checking modification indices"); FitBuoy_ok_initial = FALSE}

# STESSO PER fitAv...

######################################################
# MODIFICATION INDICES AND CORRELATED RESIDUALS (up to 2)
######################################################

# per provare
FitBuoy_ok_initial = FALSE


if (!FitBuoy_ok_initial) {
  print("Inspecting modification indices for buoyancy...")
  
  # Get modification indices (MI > 10 as threshold for inspection)
  mi_buoy <- modindices(fitBuoy, sort = TRUE, maximum.number = 10)
  print(mi_buoy)
  
  # Identify correlated residuals (theoretical candidates: similar content)
  corr_res_buoy <- mi_buoy[mi_buoy$lhs != mi_buoy$rhs & 
                             mi_buoy$op == "~~" & 
                             mi_buoy$mi > 10, ]
  
  if (nrow(corr_res_buoy) > 0 && nrow(corr_res_buoy) <= 2) {
    # specify here the justofied pairs
    # es. buoy1~~buoy2 (similar content), buoy3~~buoy4 (similar difficulty)
    justified_pairs <- c("buoy1 ~~ buoy2", "buoy3 ~~ buoy4")  # modify here
    
    if (any(corr_res_buoy$rel == justified_pairs[1]) || 
        any(corr_res_buoy$rel == justified_pairs[2])) {
      
      # Re-specify model with up to 2 correlated residuals
      respec_model_buoy <- paste0(model_buoy, "\n",
                                  justified_pairs[1],  # First justified pair
                                  ifelse(nrow(corr_res_buoy) >= 2, 
                                         paste0("\n", justified_pairs[2]), 
                                         ""))
      
      fitBuoy_respec <- cfa(respec_model_buoy, data = d_buoy, ordered = names(d_buoy))
      
      # Re-evaluate fit
      fm_respec <- fitMeasures(fitBuoy_respec, fit.measures=c("CFI","RMSEA"))
      print(paste("Respecified model fit - CFI:", round(fm_respec["cfi"], 3), 
                  "RMSEA:", round(fm_respec["rmsea"], 3)))
      
      if(fm_respec["cfi"] > 0.95 & fm_respec["rmsea"] < 0.08) {
        fitBuoy <- fitBuoy_respec  # Use respecified model
        FitBuoy_ok_initial <- TRUE
        print("Respecified model meets criteria!")
      } else {
        print("Respecified model still does not meet criteria")
      }
    } else {
      print("Suggested modifications lack theoretical justification")
    }
  } else {
    print("Too many (>2) or no suitable modification indices found")
  }
}

# Ripeti esattamente la stessa logica per avoidance...
# (fitAv_ok_initial, mi_av, corr_res_av, justified_pairs_av, etc.)

# Final fit_ok status (after possible respecification)
(fm = fitMeasures(fitBuoy, fit.measures=c("CFI","RMSEA")))
FitBuoy_ok <- isTRUE(fm["cfi"] > 0.95 && fm["rmsea"] < 0.08)

# STESSO PER fitAv_ok...

print(paste("Final buoyancy status:", ifelse(FitBuoy_ok, "OK", "EXCLUDED")))
print(paste("Final avoidance status:", ifelse(FitAv_ok, "OK", "EXCLUDED")))


# 2. Checking reliability of latent factor
(rel = semTools::reliability(fitBuoy))
if(rel["alpha.ord",] > 0.50){print("Reliability ok: proceed!"); RelBuoy_ok = T
}else {print("Stop there!")}

(rel = semTools::reliability(fitAv))
if(rel["alpha.ord",] > 0.50){print("Reliability ok: proceed!"); RelAv_ok = T
}else {print("Stop there!")}


########################################################################
# RQ2: STAGE 1 = PREDICTOR SELECTION, STAGE 2 = cs() SELECTION
########################################################################

direct_terms <- c()
if (FitBuoy_ok & RelBuoy_ok) direct_terms <- c(direct_terms, "buoyancy")
if (FitAv_ok & RelAv_ok) direct_terms <- c(direct_terms, "avoidance")

rhs_M0 <- c("worry", "sc_state", "math_ability_T0", "Gender", "(1|id)")
rhs_M1 <- c("worry", "sc_state", direct_terms, "math_ability_T0", "Gender", "(1|id)")
rhs_M2 <- c("worry", "sc_state", direct_terms, "math_anxiety", "sc_trait", "math_ability_T0", "Gender", "(1|id)")

formula_M0 <- as.formula(paste("choice_ord ~", paste(rhs_M0, collapse = " + ")))
formula_M1 <- as.formula(paste("choice_ord ~", paste(rhs_M1, collapse = " + ")))
formula_M2 <- as.formula(paste("choice_ord ~", paste(rhs_M2, collapse = " + ")))

# path a models
fit_a_worry <- brm(
  worry ~ math_anxiety + math_ability_T0 + Gender + (1 | id),
  data = d,
  family = gaussian(),
  iter = 1000,
  warmup = 500,
  chains = 4,
  cores = 4,
  seed = 220426
)

fit_a_scstate <- brm(
  sc_state ~ sc_trait + math_ability_T0 + Gender + (1 | id),
  data = d,
  family = gaussian(),
  iter = 1000,
  warmup = 500,
  chains = 4,
  cores = 4,
  seed = 220426
)

# ---------------------------
# Stage 1: M0 vs M1 vs M2
# ---------------------------

fit_M0 <- fit_choice_model(formula_M0, d)
fit_M0 <- add_criterion(fit_M0, c("loo", "waic"))

selected_stage1_model <- fit_M0
selected_stage1_name  <- "M0"
stage1_terms          <- rhs_M0
mediation_type        <- "full"

fit_M1 <- NULL
fit_M2 <- NULL

if (length(direct_terms) > 0) {
  fit_M1 <- fit_choice_model(formula_M1, d)
  fit_M1 <- add_criterion(fit_M1, c("loo", "waic"))
  
  comparison_M0_M1 <- loo_compare(fit_M0, fit_M1)
  print(comparison_M0_M1, simplify = FALSE)
  
  cmp_M0_M1 <- compare_by_elpd(comparison_M0_M1)
  
  if (cmp_M0_M1$best == "fit_M1" && cmp_M0_M1$clear) {
    selected_stage1_model <- fit_M1
    selected_stage1_name  <- "M1"
    stage1_terms          <- rhs_M1
    mediation_type        <- "full"
    
    fit_M2 <- fit_choice_model(formula_M2, d)
    fit_M2 <- add_criterion(fit_M2, c("loo", "waic"))
    
    comparison_M1_M2 <- loo_compare(fit_M1, fit_M2)
    print(comparison_M1_M2, simplify = FALSE)
    
    cmp_M1_M2 <- compare_by_elpd(comparison_M1_M2)
    
    if (cmp_M1_M2$best == "fit_M2" && cmp_M1_M2$clear) {
      selected_stage1_model <- fit_M2
      selected_stage1_name  <- "M2"
      stage1_terms          <- rhs_M2
      mediation_type        <- "partial"
    }
  }
}

# ---------------------------
# Stage 2 (category-specific effects) MODEL COMPARISON
# ---------------------------
# The code below is only an example showing how to compare:
# 1) the Stage-1 selected ordinal model
# 2) the same model with one predictor specified with cs()
#
# In the final analysis script, all relevant candidate predictors
# should be inserted one by one (or according to the pre-registered
# sequence), and all corresponding models should then be compared
# using LOO.

# starting model from Stage 1
fit_final_ord <- selected_stage1_model

# example: proportional-odds version of the selected Stage-1 model
formula_stage2_base <- as.formula(
  paste("choice_ord ~", paste(stage1_terms, collapse = " + "))
)

# example: same model, but with one predictor specified as category-specific
formula_stage2_cs_worry <- as.formula(
  paste(
    "choice_ord ~",
    paste(replace(stage1_terms, stage1_terms == "worry", "cs(worry)"), collapse = " + ")
  )
)

# fit the example category-specific model
fit_stage2_cs_worry <- fit_choice_model(formula_stage2_cs_worry, d)
fit_stage2_cs_worry <- add_criterion(fit_stage2_cs_worry, c("loo", "waic"))

# compare the Stage-1 selected model with the cs() extension
comparison_cs_example <- loo_compare(fit_final_ord, fit_stage2_cs_worry)
print(comparison_cs_example, simplify = FALSE)

cmp_cs_example <- compare_by_elpd(comparison_cs_example)

if (cmp_cs_example$best == "fit_stage2_cs_worry" && cmp_cs_example$clear) {
  final_choice_model <- fit_stage2_cs_worry
  selected_structure_label <- "cs(worry)"
} else {
  final_choice_model <- fit_final_ord
  selected_structure_label <- "proportional_odds"
}

cat("Selected Stage-1 model:", selected_stage1_name, "\n")
cat("Selected mediation type:", mediation_type, "\n")
cat("Selected ordinal structure:", selected_structure_label, "\n")

########################################################################
# INDIRECT EFFECTS (a*b)
########################################################################

draws_a_worry <- as_draws_df(fit_a_worry)
draws_a_scstate <- as_draws_df(fit_a_scstate)

a_ma_to_worry <- draws_a_worry$b_math_anxiety
a_sct_to_scstate <- draws_a_scstate$b_sc_trait

b_worry_list <- extract_univariate_choice_draws(final_choice_model, "worry")
b_scstate_list <- extract_univariate_choice_draws(final_choice_model, "sc_state")

indirect_list <- list()

for (nm in names(b_worry_list)) {
  indirect_list[[paste0("math_anxiety -> worry -> ", nm)]] <-
    summarise_draws(a_ma_to_worry * b_worry_list[[nm]], direction = "negative", label = paste0("MA_via_worry::", nm))
}

for (nm in names(b_scstate_list)) {
  indirect_list[[paste0("sc_trait -> sc_state -> ", nm)]] <-
    summarise_draws(a_sct_to_scstate * b_scstate_list[[nm]], direction = "positive", label = paste0("SC_via_sc_state::", nm))
}

results_indirect <- bind_rows(indirect_list) %>%
  mutate(
    Model = "ordinal",
    Selected_stage1_model = selected_stage1_name,
    Mediation_type = mediation_type,
    Selected_structure = selected_structure_label,
    p_one_tailed = 1 - Post.Prob,
    p_fdr = p.adjust(p_one_tailed, method = "fdr")
  )

########################################################################
# FINAL INFERENCE ON DIRECT EFFECTS
########################################################################

choice_effects <- list()

# primary state effects
for (nm in names(extract_univariate_choice_draws(final_choice_model, "worry"))) {
  choice_effects[[paste0("worry::", nm)]] <- summarise_draws(
    extract_univariate_choice_draws(final_choice_model, "worry")[[nm]],
    direction = "negative",
    label = paste0("worry::", nm)
  )
}

for (nm in names(extract_univariate_choice_draws(final_choice_model, "sc_state"))) {
  choice_effects[[paste0("sc_state::", nm)]] <- summarise_draws(
    extract_univariate_choice_draws(final_choice_model, "sc_state")[[nm]],
    direction = "positive",
    label = paste0("sc_state::", nm)
  )
}

# conditional trait predictors
if (include_buoy) {
  for (nm in names(extract_univariate_choice_draws(final_choice_model, "buoyancy"))) {
    choice_effects[[paste0("buoyancy::", nm)]] <- summarise_draws(
      extract_univariate_choice_draws(final_choice_model, "buoyancy")[[nm]],
      direction = "positive",
      label = paste0("buoyancy::", nm)
    )
  }
}

if (include_avoid) {
  for (nm in names(extract_univariate_choice_draws(final_choice_model, "avoidance"))) {
    choice_effects[[paste0("avoidance::", nm)]] <- summarise_draws(
      extract_univariate_choice_draws(final_choice_model, "avoidance")[[nm]],
      direction = "negative",
      label = paste0("avoidance::", nm)
    )
  }
}

# direct trait paths only in the partial model
if (mediation_type == "partial") {
  for (nm in names(extract_univariate_choice_draws(final_choice_model, "math_anxiety"))) {
    choice_effects[[paste0("math_anxiety::", nm)]] <- summarise_draws(
      extract_univariate_choice_draws(final_choice_model, "math_anxiety")[[nm]],
      direction = "negative",
      label = paste0("math_anxiety::", nm)
    )
  }

  for (nm in names(extract_univariate_choice_draws(final_choice_model, "sc_trait"))) {
    choice_effects[[paste0("sc_trait::", nm)]] <- summarise_draws(
      extract_univariate_choice_draws(final_choice_model, "sc_trait")[[nm]],
      direction = "positive",
      label = paste0("sc_trait::", nm)
    )
  }
}

results_rq2 <- bind_rows(choice_effects) %>%
  mutate(
    p_one_tailed = 1 - Post.Prob,
    p_fdr = p.adjust(p_one_tailed, method = "fdr"),
    Selected_stage1_model = selected_stage1_name,
    Mediation_type = mediation_type,
    Selected_structure = selected_structure_label
  )
print(results_rq2)

########################################################################
# SAVE OBJECTS
########################################################################

saveRDS(selected_stage1_model, file = file.path(Path, "selected_stage1_model.rds"))
saveRDS(final_choice_model, file = file.path(Path, "final_choice_model.rds"))
write.csv(screening_summary, file = file.path(Path, "screening_summary.csv"), row.names = FALSE)
write.csv(results_indirect, file = file.path(Path, "results_indirect.csv"), row.names = FALSE)
write.csv(results_rq2, file = file.path(Path, "results_rq2.csv"), row.names = FALSE)
