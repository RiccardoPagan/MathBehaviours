# Fit one brms model (name passed as argument), preliminary settings (2 chains, few iterations)
# Models without item-level random effects are fit on data aggregated to child x choice counts with weights():
# the likelihood is identical to the item-level one, and item-level PSIS-LOO is reconstructed in 03_summarise.R
args = commandArgs(TRUE); mod = args[1]
suppressMessages(library(brms)); source("00_prep.R")

ag = aggregate(rep(1, nrow(dm)), by = dm[, c("id", "class", "choice", zpreds, "Gender")], FUN = sum)
names(ag)[names(ag) == "x"] = "n"
ag$choice_c = factor(ag$choice, levels = c("skip", "help", "self"))
ag$choice_ord = ordered(ag$choice, levels = c("skip", "help", "self"))
ag$choice_seq = ordered(ag$choice, levels = c("self", "help", "skip"))

specs = list(
  cat       = list(f = paste("choice_c | weights(n) ~", rhs, "+ (1|id)"), fam = categorical(), dat = ag),
  cum       = list(f = paste("choice_ord | weights(n) ~", rhs, "+ (1|id)"), fam = cumulative(), dat = ag),
  acat_cs   = list(f = paste("choice_ord | weights(n) ~", rhs_cs, "+ (1|id)"), fam = acat(), dat = ag),
  seq       = list(f = paste("choice_seq | weights(n) ~", rhs, "+ (1|id)"), fam = sratio(), dat = ag),
  seq_cs    = list(f = paste("choice_seq | weights(n) ~", rhs_cs, "+ (1|id)"), fam = sratio(), dat = ag),
  cat_class = list(f = paste("choice_c | weights(n) ~", rhs, "+ (1|id) + (1|class)"), fam = categorical(), dat = ag),
  cat_item  = list(f = paste("choice_c ~", rhs, "+ (1|id) + (1|class) + (1|item)"), fam = categorical(), dat = dm)
)
s = specs[[mod]]
iter = if (mod == "cat_item") 600 else 1000
t1 = Sys.time()
fit = brm(as.formula(s$f), data = s$dat, family = s$fam, chains = 2, cores = 2, iter = iter, warmup = iter / 2,
          seed = 2026, refresh = 50, file = paste0("results/fit_", mod))
cat("\nDONE", mod, "minutes:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "\n")
print(summary(fit))
