# Fast analyses: descriptives, psychometrics, bivariate patterns, frequentist GLMMs (suppression check, clustering)
library(lme4); library(lavaan); library(semTools); library(psych); library(parallel)
source("00_prep.R")
sink("results/01_descriptives_freq.txt", split = TRUE)
options(width = 150, digits = 3)

cat("\n#### SAMPLE ####\n")
cat("Children:", length(unique(d1$id)), " classes:", length(unique(d1$class)), " schools:", length(unique(substr(d1$class, 1, 2))), "\n")
print(table(u$Gender, u$Grade))
print(table(table(d1$class)) ) # children per class (x32 rows)

cat("\n#### CHOICES ####\n")
print(round(100 * prop.table(table(d1$choice, useNA = "always")), 1))
print(round(100 * prop.table(table(d1$choice, d1$block), 2), 1))
pp = data.frame(id = sort(unique(dm$id)))
for (k in c("skip", "help", "self")) pp[[paste0("p_", k)]] = tapply(dm$choice == k, dm$id, mean)[pp$id]
pp = merge(pp, u, by = "id")
cat("All self:", sum(pp$p_self == 1), " all help:", sum(pp$p_help == 1), " all skip:", sum(pp$p_skip == 1), " never skip:", sum(pp$p_skip == 0), " never help:", sum(pp$p_help == 0), "\n")
print(summary(pp[, c("p_skip", "p_help", "p_self")]))
cat("Missing choices per child:\n"); print(table(tapply(is.na(d1$choice), d1$id, sum)))
cat("Accuracy of missing-choice rows:\n"); print(table(d1$accuracy[is.na(d1$choice)], useNA = "always"))

cat("\n#### ACCURACY BY CHOICE ####\n")
print(round(tapply(dm$accuracy, dm$choice, mean), 3))
da = subset(dm, choice != "skip"); da$help = as.numeric(da$choice == "help")
m_acc = glmer(accuracy ~ help + math_ability_T0_z + (1 | id) + (1 | item), data = da, family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("Accuracy ~ help (vs self) + ability:\n"); print(summary(m_acc)$coefficients)

cat("\n#### PREDICTOR DESCRIPTIVES (raw scales) ####\n")
print(round(psych::describe(u[, preds])[, c("n", "mean", "sd", "min", "max", "skew", "kurtosis")], 2))
cat("\nBy gender x grade (M, SD):\n")
print(aggregate(u[, preds], by = list(Gender = u$Gender, Grade = u$Grade), FUN = function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))))
cat("\nCorrelations (Pearson, child level):\n"); print(round(cor(u[, preds]), 2))

cat("\n#### PSYCHOMETRICS (item-level available only for buoyancy and avoidance) ####\n")
it = subset(d, session == 2); it = it[!duplicated(it$id), ]
bu = it[, grep("trait_buoy_", names(it))]; av = it[, grep("trait_avoidance_", names(it))]
fb = cfa(paste("buoy =~", paste(names(bu), collapse = " + ")), data = bu, ordered = names(bu))
fa = cfa(paste("av =~", paste(names(av), collapse = " + ")), data = av, ordered = names(av))
for (f in list(buoyancy = fb, avoidance = fa)) {
  print(round(fitMeasures(f, c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "rmsea.scaled", "srmr", "cfi", "rmsea")), 3))
  print(round(semTools::compRelSEM(f, ord.scale = TRUE), 3))
  print(round(standardizedSolution(f)[standardizedSolution(f)$op == "=~", "est.std"], 2))
}
cat("Ordinal alpha buoyancy / avoidance:", round(psych::alpha(psych::polychoric(bu)$rho, n.obs = nrow(bu))$total$raw_alpha, 3), "/",
    round(psych::alpha(psych::polychoric(av)$rho, n.obs = nrow(av))$total$raw_alpha, 3), "\n")

cat("\n#### BIVARIATE: child-level choice proportions x predictors ####\n")
print(round(cor(pp[, c("p_skip", "p_help", "p_self")], pp[, preds]), 2))
cat("\nVIF (child-level lm on p_self):\n")
X = pp[, zpreds]; print(round(sapply(zpreds, function(v) 1 / (1 - summary(lm(X[[v]] ~ ., data = X[, setdiff(zpreds, v)]))$r.squared)), 2))

cat("\n#### SLIDE ROBUSTNESS LM (skip=0, help=1, self=2) replicated, and p_self ####\n")
dm$score = c(skip = 0, help = 1, self = 2)[dm$choice]
pp$score = tapply(dm$score, dm$id, sum)[pp$id]
f_lm = as.formula(paste("scale(score) ~", rhs)); print(round(summary(lm(f_lm, data = pp))$coefficients, 3))
f_lm2 = as.formula(paste("scale(p_self) ~", rhs)); print(round(summary(lm(f_lm2, data = pp))$coefficients, 3))

cat("\n#### BINARY GLMMs (crossed id/class/item): full model and one-predictor-at-a-time (+ ability + gender) ####\n")
# outcomes: self vs not (stage 1 of a sequential decision), skip vs help among non-self (stage 2), skip vs not
dm$y_self = as.numeric(dm$choice == "self"); dm$y_skip = as.numeric(dm$choice == "skip")
dns = subset(dm, choice != "self"); dns$y_skip_vs_help = as.numeric(dns$choice == "skip")
specs = list()
for (y in c("y_self", "y_skip_vs_help", "y_skip")) {
  specs[[length(specs) + 1]] = list(y = y, pred = "FULL", f = paste(y, "~", rhs, "+ (1|id) + (1|class) + (1|item)"))
  specs[[length(specs) + 1]] = list(y = y, pred = "FULL_idonly", f = paste(y, "~", rhs, "+ (1|id)"))
  for (v in setdiff(zpreds, "math_ability_T0_z"))
    specs[[length(specs) + 1]] = list(y = y, pred = v, f = paste(y, "~", v, "+ math_ability_T0_z + Gender + (1|id) + (1|class) + (1|item)"))
}
cl = makeCluster(min(8, detectCores() - 2)); clusterEvalQ(cl, library(lme4)); clusterExport(cl, c("dm", "dns"))
fits = parLapply(cl, specs, function(s) {
  dat = if (s$y == "y_skip_vs_help") dns else dm
  m = glmer(as.formula(s$f), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  list(spec = s, coef = summary(m)$coefficients, vc = as.data.frame(VarCorr(m)), singular = isSingular(m))
})
stopCluster(cl)
saveRDS(fits, "results/01_glmm_fits.rds")
tab = do.call(rbind, lapply(fits, function(x) { co = x$coef; data.frame(outcome = x$spec$y, model = x$spec$pred, term = rownames(co), est = co[, 1], se = co[, 2], p = co[, 4]) }))
tab = subset(tab, term != "(Intercept)" & (model %in% c("FULL", "FULL_idonly") | term == model))
write.csv(tab, "results/01_glmm_coefs.csv", row.names = FALSE)
print(tab, row.names = FALSE)
cat("\nRandom-effect SDs (full models):\n")
for (x in fits) if (x$spec$pred %in% c("FULL", "FULL_idonly")) { cat(x$spec$y, x$spec$pred, "singular:", x$singular, "\n"); print(x$vc[, c("grp", "sdcor")]) }
sink()
