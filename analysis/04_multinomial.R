#### 04. ASSUMPTION CHECK AND MULTINOMIAL MODEL (exploratory, not preregistered) ####
# Input:  analysis/results/trials.csv (from 01_data.R)
# Output: analysis/results/04_multinomial.txt, 04_loo.csv, 04_coefs.csv, 04_coefs_class.csv, 04_acat_cs.csv
#         04_parameters.csv, 04_diagnostics.csv; analysis/figures/fig_pred.pdf/.png
# MCMC fits remain in memory only; only summaries and figures are saved.
# Paper:  Assumption Check (LOO, adjacent-category model), Multinomial Model (table + figure),
#         Robustness: class random intercept

if (basename(getwd()) == "analysis") setwd("..")
library(brms)
library(loo)
library(posterior)
library(ggplot2)

chains = 4; iter = 2000; warmup = 1000; seed = 2026
backend = "rstan"   # "cmdstanr" if CmdStan is installed (faster)

tr = read.csv("analysis/results/trials.csv")
dir.create("analysis/figures", showWarnings = FALSE)
sink("analysis/results/04_multinomial.txt", split = TRUE)
options(width = 150, digits = 3)

# counts per child x choice: same likelihood as the trial-level data, much faster (see 03)
tm = tr[!is.na(tr$choice), ]
zpreds = c("worry_z", "sc_state_z", "math_anxiety_z", "sc_trait_z", "buoyancy_z", "avoidance_z", "math_ability_T0_z")
ag = aggregate(rep(1, nrow(tm)), by = tm[, c("id", "class", "choice", zpreds, "Gender")], FUN = sum)
names(ag)[names(ag) == "x"] = "n"
ag$choice_c = factor(ag$choice, levels = c("skip", "help", "self"))
ag$choice_ord = ordered(ag$choice, levels = c("skip", "help", "self"))
rhs = paste(c(zpreds, "Gender"), collapse = " + ")
rhs_cs = paste(c(paste0("cs(", zpreds, ")"), "cs(Gender)"), collapse = " + ")

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
  write.csv(parameters, "analysis/results/04_parameters.csv", row.names = FALSE)
  write.csv(diagnostics, "analysis/results/04_diagnostics.csv", row.names = FALSE)
  print(summary(fit)); print(prior_summary(fit))
  cat("  max Rhat:", round(max(s$rhat, na.rm = TRUE), 3), " min bulk ESS:", round(min(s$ess_bulk, na.rm = TRUE)),
      " divergences:", sum(np$Value[np$Parameter == "divergent__"]), "\n")
}

#### models ####
cat("cat\n")
fit_cat = brm(as.formula(paste("choice_c | weights(n) ~", rhs, "+ (1 | id)")), data = ag, family = categorical(),
              chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_cat, "04_cat")
cat("cum\n")
fit_cum = brm(as.formula(paste("choice_ord | weights(n) ~", rhs, "+ (1 | id)")), data = ag, family = cumulative(),
              chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_cum, "04_cum")
# adjacent-category model with category-specific effects for all predictors:
# shows whether a predictor moves children towards the middle category (help)
cat("acat_cs\n")
fit_acat = brm(as.formula(paste("choice_ord | weights(n) ~", rhs_cs, "+ (1 | id)")), data = ag, family = acat(),
               chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_acat, "04_acat_cs")
cat("cat_class\n")
fit_class = brm(as.formula(paste("choice_c | weights(n) ~", rhs, "+ (1 | id) + (1 | class)")), data = ag, family = categorical(),
                chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, backend = backend)
check(fit_class, "04_cat_class")

#### assumption check: LOO ####
cat("\n#### LOO (trial level) ####\n")
pw = list(cat = loo_trials(fit_cat), cum = loo_trials(fit_cum), acat_cs = loo_trials(fit_acat))
looTab = data.frame(model = names(pw), elpd = sapply(pw, sum), se = sapply(pw, function(x) sqrt(length(x) * var(x))))
looTab$elpd_diff_vs_cat = sapply(pw, function(x) sum(x - pw$cat))
looTab$se_diff = sapply(pw, function(x) sqrt(length(x) * var(x - pw$cat)))
looTab$diff_in_se = looTab$elpd_diff_vs_cat / looTab$se_diff
print(looTab, row.names = FALSE)
write.csv(looTab, "analysis/results/04_loo.csv", row.names = FALSE)

cat("\n#### ADJACENT-CATEGORY MODEL, category-specific effects ([1] = help vs skip, [2] = self vs help) ####\n")
fe = fixef(fit_acat)
print(round(fe, 3))
write.csv(round(fe, 3), "analysis/results/04_acat_cs.csv")

#### multinomial coefficients: all three contrasts ####
# reference = skip; self vs help derived from the draws
coefs = function(fit, label) {
  dr = as_draws_df(fit)
  out = data.frame()
  for (v in c(zpreds, "GenderM")) {
    h = dr[[paste0("b_muhelp_", v)]]
    s = dr[[paste0("b_muself_", v)]]
    for (k in c("help_vs_skip", "self_vs_skip", "self_vs_help")) {
      x = switch(k, help_vs_skip = h, self_vs_skip = s, self_vs_help = s - h)
      out = rbind(out, data.frame(model = label, contrast = k, term = v, est = mean(x), lo = quantile(x, .025),
                                  hi = quantile(x, .975), pp_pos = mean(x > 0), row.names = NULL))
    }
  }
  out
}
cat("\n#### MULTINOMIAL MODEL (reference = skip; pp_pos = posterior probability of a positive effect) ####\n")
co = coefs(fit_cat, "cat")
print(co, row.names = FALSE)
write.csv(co, "analysis/results/04_coefs.csv", row.names = FALSE)
print(VarCorr(fit_cat))

cat("\n#### WITH CLASS RANDOM INTERCEPT ####\n")
co2 = coefs(fit_class, "cat_class")
co2$change_vs_cat = co2$est - co$est
print(co2, row.names = FALSE)
cat("Largest change:", round(max(abs(co2$change_vs_cat)), 3), "\n")
write.csv(co2, "analysis/results/04_coefs_class.csv", row.names = FALSE)
print(VarCorr(fit_class))

#### figure: predicted probabilities ####
# population-level predictions (child intercept = 0), other predictors at 0 (their mean),
# averaged over gender with the sample proportions
show = c(math_ability_T0_z = "Math ability (T0)", avoidance_z = "Math avoidance",
         buoyancy_z = "Academic buoyancy", worry_z = "Worry (state)")
pF = mean(tr$Gender[!duplicated(tr$id)] == "F")
grid = seq(-2, 2, by = 0.1)
pl = data.frame()
for (v in names(show)) {
  nd = data.frame(matrix(0, nrow = length(grid), ncol = length(zpreds), dimnames = list(NULL, zpreds)))
  nd[[v]] = grid
  nd$n = 1
  ndF = nd; ndF$Gender = "F"
  ndM = nd; ndM$Gender = "M"
  ep = pF * posterior_epred(fit_cat, newdata = ndF, re_formula = NA) +
    (1 - pF) * posterior_epred(fit_cat, newdata = ndM, re_formula = NA)   # draws x grid x category
  for (k in dimnames(ep)[[3]]) {
    pl = rbind(pl, data.frame(predictor = show[v], z = grid, choice = k, p = apply(ep[, , k], 2, mean),
                              lo = apply(ep[, , k], 2, quantile, .025), hi = apply(ep[, , k], 2, quantile, .975)))
  }
}
pl$predictor = factor(pl$predictor, levels = show)
pl$choice = factor(pl$choice, levels = c("skip", "help", "self"))
write.csv(pl, "analysis/results/04_pred.csv", row.names = FALSE)

g = ggplot(pl, aes(x = z, y = p, colour = choice, fill = choice)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ predictor, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Predictor (SD units)", y = "Predicted probability", colour = "Choice", fill = "Choice") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("analysis/figures/fig_pred.pdf", g, width = 9, height = 3.5)
ggsave("analysis/figures/fig_pred.png", g, width = 9, height = 3.5, dpi = 300)

sink()
