#### 02. DESCRIPTIVES, PSYCHOMETRIC SCREENING, BIVARIATE PATTERNS ####
# Input:  analysis/results/trials.csv, children.csv (from 01_data.R)
# Output: analysis/results/02_descriptives.txt (full printout) + tables tab_*.csv
# Paper:  Descriptive Statistics, Psychometric Screening, Assumption Check (bivariate table)

if (basename(getwd()) == "analysis") setwd("..")
library(lme4)
library(lavaan)
library(semTools)
library(psych)

tr = read.csv("analysis/results/trials.csv")
ch = read.csv("analysis/results/children.csv")
preds = c("worry", "sc_state", "math_anxiety", "sc_trait", "buoyancy", "avoidance", "math_ability_T0")
sink("analysis/results/02_descriptives.txt", split = TRUE)
options(width = 150, digits = 3)
report = list()

#### sample ####
cat("\n#### SAMPLE ####\n")
cat("Children:", nrow(ch), " classes:", length(unique(ch$class)), " schools:", length(unique(ch$school)), "\n")
print(table(ch$Gender, ch$Grade))
cat("Children per class:\n"); print(summary(as.numeric(table(ch$class))))
report$sample = list(n = nrow(ch), girls = sum(ch$Gender == "F"), classes = length(unique(ch$class)),
                     schools = length(unique(ch$school)), class_range = range(table(ch$class)))

#### predictors ####
cat("\n#### PREDICTORS (raw scores) ####\n")
desc = psych::describe(ch[, preds])[, c("mean", "sd", "min", "max", "skew", "kurtosis")]
print(round(desc, 2))
write.csv(round(desc, 2), "analysis/results/tab_descriptives.csv")

bygg = aggregate(ch[, preds], by = list(Gender = ch$Gender, Grade = ch$Grade),
                 FUN = function(x) sprintf("%.2f (%.2f)", mean(x), sd(x)))
bygg$n = aggregate(ch$id, by = list(ch$Gender, ch$Grade), FUN = length)$x
print(bygg)
write.csv(bygg, "analysis/results/tab_gender_grade.csv", row.names = FALSE)

r = cor(ch[, preds])
print(round(r, 2))
write.csv(round(r, 2), "analysis/results/tab_correlations.csv")

#### choices ####
cat("\n#### CHOICES ####\n")
print(round(100 * prop.table(table(tr$choice, useNA = "always")), 1))
cat("By block (% of valid choices):\n"); print(round(100 * prop.table(table(tr$choice, tr$block), 2), 1))
cat("All self:", sum(ch$p_self == 1), " all help:", sum(ch$p_help == 1), " all skip:", sum(ch$p_skip == 1),
    " never skip:", sum(ch$n_skip == 0), " never help:", sum(ch$n_help == 0), "\n")
print(summary(ch[, c("p_skip", "p_help", "p_self")]))
cat("Missing choices per child:\n"); print(table(ch$n_missing))
cat("Accuracy of trials with missing choice:\n"); print(table(tr$accuracy[is.na(tr$choice)], useNA = "always"))
report$choices = list(percent = 100 * prop.table(table(tr$choice, useNA = "always")),
                     block = 100 * prop.table(table(tr$choice, tr$block), 2),
                     always = c(self = sum(ch$p_self == 1), help = sum(ch$p_help == 1), skip = sum(ch$p_skip == 1)),
                     never = c(skip = sum(ch$n_skip == 0), help = sum(ch$n_help == 0)),
                     median_skip = median(ch$p_skip), missing = sum(is.na(tr$choice)),
                     no_missing = sum(ch$n_missing == 0), max_missing = max(ch$n_missing))

#### accuracy by choice ####
cat("\n#### ACCURACY BY CHOICE ####\n")
tm = tr[!is.na(tr$choice), ]
print(round(tapply(tm$accuracy, tm$choice, mean), 3))
# does the hint help, at equal ability? (trials with help or self only)
da = tm[tm$choice != "skip", ]
da$help = as.numeric(da$choice == "help")
fit_acc = glmer(accuracy ~ help + math_ability_T0_z + (1 | id) + (1 | item), data = da, family = binomial,
                control = glmerControl(optimizer = "bobyqa"))
print(summary(fit_acc)$coefficients)
report$accuracy = list(by_choice = tapply(tm$accuracy, tm$choice, mean), coefficients = summary(fit_acc)$coefficients)

#### psychometric screening (preregistered: buoyancy and avoidance only) ####
cat("\n#### PSYCHOMETRIC SCREENING ####\n")
# Preregistered criteria: CFI > .95 and RMSEA < .08. Both the scaled (robust, WLSMV) and the
# unscaled indices are printed: the paper must state which one the criterion refers to.
bu = ch[, paste0("trait_buoy_", 1:4)]
av = ch[, paste0("trait_avoidance_", 1:10)]
fit_bu = cfa(paste("buoy =~", paste(names(bu), collapse = " + ")), data = bu, ordered = names(bu))
fit_av = cfa(paste("av =~", paste(names(av), collapse = " + ")), data = av, ordered = names(av))

idx = c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "rmsea.scaled", "srmr", "cfi", "rmsea")
report$psych = list()
for (v in c("buoyancy", "avoidance")) {
  f = if (v == "buoyancy") fit_bu else fit_av
  print(round(fitMeasures(f, idx), 3))
  std = standardizedSolution(f)
  cat("Standardised loadings:", round(std$est.std[std$op == "=~"], 2), "\n")
  cat("Omega (ordinal):", round(compRelSEM(f, ord.scale = TRUE), 3), "\n")
  report$psych[[v]] = list(indices = fitMeasures(f, idx), loadings = std$est.std[std$op == "=~"],
                           omega = as.numeric(compRelSEM(f, ord.scale = TRUE)))
}
cat("Ordinal alpha buoyancy:", round(psych::alpha(polychoric(bu)$rho, n.obs = nrow(bu))$total$raw_alpha, 3), "\n")
cat("Ordinal alpha avoidance:", round(psych::alpha(polychoric(av)$rho, n.obs = nrow(av))$total$raw_alpha, 3), "\n")
cat("Buoyancy-avoidance correlation (sum scores):", round(cor(ch$buoyancy, ch$avoidance), 3), "\n")
report$psych$buoyancy$alpha = psych::alpha(polychoric(bu)$rho, n.obs = nrow(bu))$total$raw_alpha
report$psych$avoidance$alpha = psych::alpha(polychoric(av)$rho, n.obs = nrow(av))$total$raw_alpha
report$psych$r = cor(ch$buoyancy, ch$avoidance)

# If the avoidance model fails the criterion: at most two correlated residuals, only between items
# with clear overlap in content or wording (preregistered). Inspect the largest modification indices
# here, check the item wording, and only then set the residuals in 03_preregistered.R.
cat("\nLargest modification indices (avoidance):\n")
mi = modindices(fit_av, sort. = TRUE)
print(head(mi[mi$op == "~~", c("lhs", "op", "rhs", "mi", "epc")], 10))

#### bivariate: per-child proportions x predictors ####
cat("\n#### BIVARIATE: per-child choice proportions x predictors ####\n")
rb = cor(ch[, c("p_skip", "p_help", "p_self")], ch[, preds])
print(round(rb, 2))
write.csv(round(rb, 2), "analysis/results/tab_bivariate.csv")

cat("\nVariance inflation factors (standardised predictors):\n")
X = ch[, paste0(preds, "_z")]
vif = c()
for (v in names(X)) vif[v] = 1 / (1 - summary(lm(X[[v]] ~ ., data = X[, names(X) != v]))$r.squared)
print(round(vif, 2))
report$vif = vif
saveRDS(report, "analysis/results/02_report.rds")  # numeric summaries only, no fitted models

sink()
