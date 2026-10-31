#### 03. PREREGISTERED PIPELINE (ordinal) ####
# Input:  analysis/results/trials.csv, children.csv (from 01_data.R)
# Output: analysis/results/03_preregistered.txt, 03_selection.csv, 03_effects.csv, 03_indirect.csv
#         03_parameters.csv, 03_diagnostics.csv; MCMC fits remain in memory only
# Paper:  Ordinal Models (Table "Final Ordinal Model")
#
# Steps, as in the preregistration:
#   1. Psychometric screening of buoyancy and avoidance (CFA; |r| < .70)
#   2. Stage 1: M0 (state + covariates) -> M1 (+ buoyancy, avoidance) -> M2 (+ math anxiety, self-competence)
#   3. Stage 2: category-specific effects, one predictor at a time, each vs its predecessor
#   4. Final model: one-tailed p (from posterior mean and SD) + FDR; indirect effects a*b
# A model is preferred only if ELPD_diff > 2 SE (PSIS-LOO at the trial level); otherwise selection stops.

if (basename(getwd()) == "analysis") setwd("..")
library(brms)
library(lavaan)
library(loo)
library(posterior)

#### DECISIONS TO SET BEFORE RUNNING ####
# (a) which CFI/RMSEA the preregistered criterion refers to: scaled (robust) or unscaled.
#     With the scaled indices avoidance fails (CFI = .92); with the unscaled ones it passes (CFI = .97).
fit_index = c("cfi.scaled", "rmsea.scaled")   # or c("cfi", "rmsea")
# (b) if avoidance fails: at most two correlated residuals, ONLY between items with clear overlap
#     in content or wording (check the wording; candidates in 02_descriptives.txt). Leave empty otherwise.
av_residuals = c()   # e.g. c("trait_avoidance_5 ~~ trait_avoidance_6")
# (c) order in which category-specific effects are tried in Stage 2 (primary predictors only;
#     those not in the selected model are skipped). The preregistration does not fix the order.
cs_order = c("worry_z", "sc_state_z", "buoyancy_z", "avoidance_z")

# MCMC settings (brms default priors: flat on regression coefficients)
chains = 4; iter = 2000; warmup = 1000; seed = 2026
backend = "rstan"   # "cmdstanr" if CmdStan is installed (faster)

tr = read.csv("analysis/results/trials.csv")
ch = read.csv("analysis/results/children.csv")
sink("analysis/results/03_preregistered.txt", split = TRUE)
options(width = 150, digits = 3)

#### 1. psychometric screening ####
bu = ch[, paste0("trait_buoy_", 1:4)]
av = ch[, paste0("trait_avoidance_", 1:10)]
fit_bu = cfa(paste("buoy =~", paste(names(bu), collapse = " + ")), data = bu, ordered = names(bu))
fit_av = cfa(paste(c(paste("av =~", paste(names(av), collapse = " + ")), av_residuals), collapse = "\n"),
             data = av, ordered = names(av))
fm_bu = fitMeasures(fit_bu, fit_index)
fm_av = fitMeasures(fit_av, fit_index)
include_buoy = fm_bu[1] > .95 & fm_bu[2] < .08
include_avoid = fm_av[1] > .95 & fm_av[2] < .08
r_ba = cor(ch$buoyancy, ch$avoidance)
cat("Buoyancy CFA:", round(fm_bu, 3), "-> included:", include_buoy, "\n")
cat("Avoidance CFA:", round(fm_av, 3), " residuals:", av_residuals, "-> included:", include_avoid, "\n")
cat("Buoyancy-avoidance r =", round(r_ba, 3), "\n")
if (abs(r_ba) > .70) stop("|r| > .70: the preregistration requires separate models for buoyancy and avoidance")

#### data: counts per child x choice ####
# The likelihood of the child x choice counts (weights = number of trials) is identical to the
# trial-level likelihood, but sampling is ~15 times faster. Trial-level LOO is reconstructed below.
tm = tr[!is.na(tr$choice), ]
zpreds = c("worry_z", "sc_state_z", "math_anxiety_z", "sc_trait_z", "buoyancy_z", "avoidance_z", "math_ability_T0_z")
ag = aggregate(rep(1, nrow(tm)), by = tm[, c("id", "choice", zpreds, "Gender")], FUN = sum)
names(ag)[names(ag) == "x"] = "n"
ag$choice_ord = ordered(ag$choice, levels = c("skip", "help", "self"))
cat("Trials analysed:", sum(ag$n), " children:", length(unique(ag$id)), "\n")

fit_ord = function(terms, name) {
  f = as.formula(paste("choice_ord | weights(n) ~", paste(terms, collapse = " + "), "+ (1 | id)"))
  brm(f, data = ag, family = cumulative("logit"), chains = chains, iter = iter, warmup = warmup,
      cores = chains, seed = seed, backend = backend)
}

# trial-level PSIS-LOO from a fit on counts: each row's log-lik is n * log p, so divide by n
# and repeat the pointwise elpd n times (the n trials of a row are identical)
loo_trials = function(fit) {
  ll = log_lik(fit)
  n = fit$data$n
  ll = sweep(ll, 2, n, "/")
  lo = loo(ll, r_eff = relative_eff(exp(ll), chain_id = rep(1:chains, each = iter - warmup)))
  cat("  Pareto k > 0.7 (trials):", sum(rep(lo$diagnostics$pareto_k, n) > 0.7), "\n")
  rep(lo$pointwise[, "elpd_loo"], n)
}

parameters = diagnostics = data.frame()
check = function(fit, name) {
  dr = as_draws_df(fit, variable = "^(b|sd_|cor_|sigma|disc)", regex = TRUE)
  s = as.data.frame(summarise_draws(dr, "mean", "sd", ~quantile(.x, c(.025, .975)), "rhat", "ess_bulk", "ess_tail"))
  np = nuts_params(fit)
  parameters <<- rbind(parameters, data.frame(model = name, s))
  diagnostics <<- rbind(diagnostics, data.frame(model = name, draws = ndraws(fit), chains = nchains(fit),
                     max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
                     min_ess_tail = min(s$ess_tail, na.rm = TRUE),
                     divergences = sum(np$Value[np$Parameter == "divergent__"])))
  write.csv(parameters, "analysis/results/03_parameters.csv", row.names = FALSE)
  write.csv(diagnostics, "analysis/results/03_diagnostics.csv", row.names = FALSE)
  print(summary(fit)); print(prior_summary(fit))
  cat("  max Rhat:", round(max(s$rhat, na.rm = TRUE), 3), " min bulk ESS:", round(min(s$ess_bulk, na.rm = TRUE)),
      " divergences:", sum(np$Value[np$Parameter == "divergent__"]), "\n")
}

selection = data.frame()
compare = function(pw_simple, pw_complex, step) {
  dd = pw_complex - pw_simple
  out = data.frame(step = step, elpd_diff = sum(dd), se_diff = sqrt(length(dd) * var(dd)))
  out$clear = out$elpd_diff > 2 * out$se_diff
  print(out)
  out
}

#### 2. Stage 1: M0 -> M1 -> M2 ####
cov = c("math_ability_T0_z", "Gender")
direct = c("buoyancy_z", "avoidance_z")[c(include_buoy, include_avoid)]
terms_M0 = c("worry_z", "sc_state_z", cov)
terms_M1 = c("worry_z", "sc_state_z", direct, cov)
terms_M2 = c("worry_z", "sc_state_z", direct, "math_anxiety_z", "sc_trait_z", cov)

cat("\n#### STAGE 1 ####\nM0\n")
fit_M0 = fit_ord(terms_M0, "M0"); check(fit_M0, "03_M0"); pw_M0 = loo_trials(fit_M0)
selected = "M0"; terms_sel = terms_M0; fit_sel = fit_M0; pw_sel = pw_M0

if (length(direct) > 0) {
  cat("M1\n")
  fit_M1 = fit_ord(terms_M1, "M1"); check(fit_M1, "03_M1"); pw_M1 = loo_trials(fit_M1)
  cmp = compare(pw_M0, pw_M1, "M1 vs M0"); selection = rbind(selection, cmp)
  if (cmp$clear) {
    selected = "M1"; terms_sel = terms_M1; fit_sel = fit_M1; pw_sel = pw_M1
    cat("M2\n")
    fit_M2 = fit_ord(terms_M2, "M2"); check(fit_M2, "03_M2"); pw_M2 = loo_trials(fit_M2)
    cmp = compare(pw_M1, pw_M2, "M2 vs M1"); selection = rbind(selection, cmp)
    if (cmp$clear) { selected = "M2"; terms_sel = terms_M2; fit_sel = fit_M2; pw_sel = pw_M2 }
  }
}
cat("Selected Stage-1 model:", selected, "\n")

#### 3. Stage 2: category-specific effects, one at a time ####
cat("\n#### STAGE 2 ####\n")
cs_done = c()
for (v in intersect(cs_order, terms_sel)) {
  terms_new = terms_sel
  terms_new[terms_new == v] = paste0("cs(", v, ")")
  name = paste0(selected, "_cs_", paste(c(cs_done, v), collapse = "_"))
  cat(name, "\n")
  fit_new = fit_ord(terms_new, name); check(fit_new, paste0("03_", name)); pw_new = loo_trials(fit_new)
  cmp = compare(pw_sel, pw_new, paste0("cs(", v, ") vs previous")); selection = rbind(selection, cmp)
  if (!cmp$clear) break
  terms_sel = terms_new; fit_sel = fit_new; pw_sel = pw_new; cs_done = c(cs_done, v)
}
cat("Final model terms:", terms_sel, "\n")
write.csv(selection, "analysis/results/03_selection.csv", row.names = FALSE)
saveRDS(list(fit_index = fit_index, av_residuals = av_residuals, include_buoyancy = include_buoy,
             include_avoidance = include_avoid, stage1 = selected, terms = terms_sel,
             final_model = paste0("03_", selected, if (length(cs_done)) paste0("_cs_", paste(cs_done, collapse = "_")) else ""),
             chains = chains, iter = iter, warmup = warmup, seed = seed, backend = backend,
             R = as.character(getRversion()), brms = as.character(packageVersion("brms"))),
        "analysis/results/03_spec.rds")

#### 4. final model: directional tests ####
# preregistered "pragmatic frequentist" inference: one-tailed p from posterior mean and SD, FDR
# across the primary effects retained. Positive = more engaged choice (towards self).
cat("\n#### FINAL MODEL ####\n")
print(summary(fit_sel))
expected = c(worry_z = -1, sc_state_z = 1, buoyancy_z = 1, avoidance_z = -1, math_anxiety_z = -1, sc_trait_z = 1)
dr = as_draws_df(fit_sel)
eff = data.frame()
for (v in intersect(names(expected), gsub("cs\\(|\\)", "", terms_sel))) {
  par = grep(paste0("^b(cs)?_", v, "(\\[|$)"), names(dr), value = TRUE)
  for (p in par) {
    x = dr[[p]]
    z = mean(x) / sd(x)
    eff = rbind(eff, data.frame(parameter = p, est = mean(x), se = sd(x), lo = quantile(x, .025), hi = quantile(x, .975),
                                p_one = ifelse(expected[v] < 0, pnorm(z), 1 - pnorm(z)), row.names = NULL))
  }
}
eff$p_fdr = p.adjust(eff$p_one, method = "fdr")
print(eff)
write.csv(eff, "analysis/results/03_effects.csv", row.names = FALSE)

#### indirect effects: trait -> state (path a) x state -> choice (path b) ####
# path a: one row per child (state measures do not vary across trials). Fitting it separately from
# the choice model is equivalent to the joint model with set_rescor(FALSE): no shared parameters.
fit_a_worry = brm(worry_z ~ math_anxiety_z + math_ability_T0_z + Gender, data = ch, chains = chains, iter = iter,
                  warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_a_worry, "03_a_worry")
fit_a_scstate = brm(sc_state_z ~ sc_trait_z + math_ability_T0_z + Gender, data = ch, chains = chains, iter = iter,
                    warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_a_scstate, "03_a_scstate")
a_worry = as_draws_df(fit_a_worry)$b_math_anxiety_z
a_scstate = as_draws_df(fit_a_scstate)$b_sc_trait_z

ind = data.frame()
for (p in grep("^b(cs)?_(worry_z|sc_state_z)(\\[|$)", names(dr), value = TRUE)) {
  a = if (grepl("worry", p)) a_worry else a_scstate
  x = a * dr[[p]]
  z = mean(x) / sd(x)
  dir = if (grepl("worry", p)) -1 else 1
  ind = rbind(ind, data.frame(path = paste(ifelse(dir < 0, "math_anxiety -> worry ->", "sc_trait -> sc_state ->"), p),
                              est = mean(x), se = sd(x), lo = quantile(x, .025), hi = quantile(x, .975),
                              p_one = ifelse(dir < 0, pnorm(z), 1 - pnorm(z)), row.names = NULL))
}
ind$p_fdr = p.adjust(ind$p_one, method = "fdr")
cat("Path a (math anxiety -> worry):", round(mean(a_worry), 3), " (self-competence -> perceived competence):", round(mean(a_scstate), 3), "\n")
print(ind)
write.csv(ind, "analysis/results/03_indirect.csv", row.names = FALSE)

sink()
