#### 05. ROBUSTNESS (frequentist) ####
# Input:  analysis/results/trials.csv, children.csv (from 01_data.R)
# Output: analysis/results/05_robustness.txt, 05_lm.csv, 05_glmm.csv
# Paper:  Robustness: Linear Model, Class and Item Random Intercepts, Skip vs Help, Buoyancy

if (basename(getwd()) == "analysis") setwd("..")
library(lme4)
library(parallel)

tr = read.csv("analysis/results/trials.csv")
ch = read.csv("analysis/results/children.csv")
zpreds = c("worry_z", "sc_state_z", "math_anxiety_z", "sc_trait_z", "buoyancy_z", "avoidance_z", "math_ability_T0_z")
rhs = paste(c(zpreds, "Gender"), collapse = " + ")
sink("analysis/results/05_robustness.txt", split = TRUE)
options(width = 150, digits = 3)

#### linear models on per-child scores ####
# score: skip = 0, help = 1, self = 2, summed over valid trials (assumes the order skip < help < self);
# proportion of self choices: assumes no order between skip and help
ch$score = ch$n_help + 2 * ch$n_self
fit_score = lm(as.formula(paste("scale(score) ~", rhs)), data = ch)
fit_pself = lm(as.formula(paste("scale(p_self) ~", rhs)), data = ch)
cat("\n#### LM: choice score ####\n"); print(summary(fit_score)); print(round(confint(fit_score), 2))
cat("\n#### LM: proportion of self choices ####\n"); print(summary(fit_pself)); print(round(confint(fit_pself), 2))
lmtab = rbind(data.frame(outcome = "score", term = names(coef(fit_score)), summary(fit_score)$coefficients, confint(fit_score), row.names = NULL),
              data.frame(outcome = "p_self", term = names(coef(fit_pself)), summary(fit_pself)$coefficients, confint(fit_pself), row.names = NULL))
names(lmtab)[3:8] = c("est", "se", "t", "p", "lo", "hi")
write.csv(lmtab, "analysis/results/05_lm.csv", row.names = FALSE)
lmfit = data.frame()
for (y in c("score", "p_self")) {
  s = summary(if (y == "score") fit_score else fit_pself)
  f = s$fstatistic
  lmfit = rbind(lmfit, data.frame(outcome = y, r2 = s$r.squared, adj_r2 = s$adj.r.squared,
                 F = unname(f[1]), df1 = unname(f[2]), df2 = unname(f[3]), p = pf(f[1], f[2], f[3], lower.tail = FALSE)))
}
write.csv(lmfit, "analysis/results/05_lm_fit.csv", row.names = FALSE)

#### binary mixed logistic models ####
# self vs not self; skip vs help (trials without self choice); skip vs not skip.
# For each: all predictors with child, class and item intercepts ("FULL"), the same with the child
# intercept only ("FULL_idonly"), and one attitude at a time (+ ability + gender).
tm = tr[!is.na(tr$choice), ]
tm$y_self = as.numeric(tm$choice == "self")
tm$y_skip = as.numeric(tm$choice == "skip")
tns = tm[tm$choice != "self", ]
tns$y_skip_vs_help = as.numeric(tns$choice == "skip")

specs = list()
for (y in c("y_self", "y_skip_vs_help", "y_skip")) {
  specs[[length(specs) + 1]] = list(y = y, model = "FULL", f = paste(y, "~", rhs, "+ (1|id) + (1|class) + (1|item)"))
  specs[[length(specs) + 1]] = list(y = y, model = "FULL_idonly", f = paste(y, "~", rhs, "+ (1|id)"))
  for (v in setdiff(zpreds, "math_ability_T0_z"))
    specs[[length(specs) + 1]] = list(y = y, model = v, f = paste(y, "~", v, "+ math_ability_T0_z + Gender + (1|id) + (1|class) + (1|item)"))
}

cl = makeCluster(min(8, detectCores() - 2))
clusterEvalQ(cl, library(lme4))
clusterExport(cl, c("tm", "tns"))
fits = parLapply(cl, specs, function(s) {
  dat = if (s$y == "y_skip_vs_help") tns else tm
  m = glmer(as.formula(s$f), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  list(spec = s, coef = summary(m)$coefficients, vc = as.data.frame(VarCorr(m)), singular = isSingular(m),
       warn = m@optinfo$conv$lme4$messages)
})
stopCluster(cl)
random = convergence = data.frame()
for (x in fits) {
  convergence = rbind(convergence, data.frame(outcome = x$spec$y, model = x$spec$model,
                       singular = x$singular, warnings = paste(x$warn, collapse = "; ")))
  if (x$spec$model %in% c("FULL", "FULL_idonly"))
    random = rbind(random, data.frame(outcome = x$spec$y, model = x$spec$model,
                                     group = x$vc$grp, sd = x$vc$sdcor))
}
write.csv(random, "analysis/results/05_random_effects.csv", row.names = FALSE)
write.csv(convergence, "analysis/results/05_convergence.csv", row.names = FALSE)

tab = data.frame()
for (x in fits) {
  co = x$coef
  keep = rownames(co) != "(Intercept)" & (x$spec$model %in% c("FULL", "FULL_idonly") | rownames(co) == x$spec$model)
  tab = rbind(tab, data.frame(outcome = x$spec$y, model = x$spec$model, term = rownames(co)[keep],
                              est = co[keep, 1], se = co[keep, 2], p = co[keep, 4], row.names = NULL))
}
cat("\n#### BINARY GLMMs ####\n")
print(tab, row.names = FALSE)
write.csv(tab, "analysis/results/05_glmm.csv", row.names = FALSE)

cat("\nRandom-effect SDs and convergence (full models):\n")
for (x in fits) if (x$spec$model %in% c("FULL", "FULL_idonly")) {
  cat(x$spec$y, x$spec$model, " singular:", x$singular, " warnings:", x$warn, "\n")
  print(x$vc[, c("grp", "sdcor")])
}
cat("\nAny convergence warning in the other models:", any(sapply(fits, function(x) length(x$warn) > 0)), "\n")

sink()
