# Summarise preliminary brms fits: item-level LOO, diagnostics, coefficients (all contrasts), average marginal effects
suppressMessages({library(brms); library(loo); library(posterior)})
source("00_prep.R")
sink("results/03_summary.txt", split = TRUE); options(width = 150, digits = 3)
mods = c("cat", "cum", "acat_cs", "seq", "seq_cs", "cat_class", "cat_item")
fits = list(); for (m in mods) { f = paste0("results/fit_", m, ".rds"); if (file.exists(f)) fits[[m]] = readRDS(f) }
cat("Available fits:", names(fits), "\n")

cat("\n#### DIAGNOSTICS ####\n")
diag = do.call(rbind, lapply(names(fits), function(m) {
  fit = fits[[m]]; sm = summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  np = nuts_params(fit)
  data.frame(model = m, max_rhat = max(sm$rhat, na.rm = T), min_ess_bulk = min(sm$ess_bulk, na.rm = T),
             min_ess_tail = min(sm$ess_tail, na.rm = T), divergences = sum(np$Value[np$Parameter == "divergent__"]),
             draws = ndraws(fit))
}))
print(diag); write.csv(diag, "results/03_diagnostics.csv", row.names = FALSE)

cat("\n#### ITEM-LEVEL PSIS-LOO ####\n")
# aggregated fits: log_lik returns n * log p for each child x choice row; per-item log p = ll / n,
# PSIS is identical for the n identical items, so item-level pointwise elpd = rep(elpd_row, n)
pw = list(); kbad = c()
for (m in names(fits)) {
  fit = fits[[m]]; ll = log_lik(fit)
  n = if (m == "cat_item") rep(1, ncol(ll)) else fit$data$n
  ll = sweep(ll, 2, n, "/")
  lo = loo(ll, r_eff = relative_eff(exp(ll), chain_id = rep(1:nchains(fit), each = niterations(fit))))
  pw[[m]] = rep(lo$pointwise[, "elpd_loo"], n); kbad[m] = sum(rep(lo$diagnostics$pareto_k, n) > 0.7)
}
N = sapply(pw, length); print(N)
best = names(which.max(sapply(pw, sum)))
looTab = data.frame(model = names(pw), elpd_loo = sapply(pw, sum), se = sapply(pw, function(x) sqrt(length(x) * var(x))),
                    elpd_diff_vs_best = sapply(pw, function(x) sum(x - pw[[best]])),
                    se_diff = sapply(pw, function(x) sqrt(length(x) * var(x - pw[[best]]))), n_k_gt_0.7 = kbad)
looTab$diff_in_SE = looTab$elpd_diff_vs_best / looTab$se_diff
looTab = looTab[order(-looTab$elpd_loo), ]; print(looTab, row.names = FALSE)
write.csv(looTab, "results/03_loo.csv", row.names = FALSE)
cat("Note: cat_item includes an item random effect (crossed), so it is not a pure 'response structure' comparison.\n")

summ = function(x) c(mean = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975)), pp_pos = mean(x > 0))

cat("\n#### MULTINOMIAL COEFFICIENTS (logit, reference = skip; self vs help derived) ####\n")
for (m in intersect(c("cat", "cat_class", "cat_item"), names(fits))) {
  dr = as_draws_df(fits[[m]]); terms = c(zpreds, "GenderM")
  tab = do.call(rbind, lapply(terms, function(v) {
    h = dr[[paste0("b_muhelp_", v)]]; s = dr[[paste0("b_muself_", v)]]
    rbind(data.frame(model = m, contrast = "help_vs_skip", term = v, t(summ(h))),
          data.frame(model = m, contrast = "self_vs_skip", term = v, t(summ(s))),
          data.frame(model = m, contrast = "self_vs_help", term = v, t(summ(s - h))))
  }))
  tab = tab[order(tab$contrast), ]; print(tab, row.names = FALSE)
  write.csv(tab, paste0("results/03_coefs_", m, ".csv"), row.names = FALSE)
}

cat("\n#### SEQUENTIAL MODEL WITH CATEGORY-SPECIFIC EFFECTS (sratio; order self < help < skip) ####\n")
cat("Stage [1] = continuing beyond 'self' (i.e. NOT self) vs self; stage [2] = skip vs help among non-self. Positive = more likely to continue.\n")
if ("seq_cs" %in% names(fits)) { fe = fixef(fits[["seq_cs"]]); print(round(fe, 3)); write.csv(fe, "results/03_coefs_seq_cs.csv") }

cat("\n#### AVERAGE MARGINAL EFFECTS on choice probabilities (per 1 SD, from -0.5 to +0.5 SD; Gender: M - F) ####\n")
nd0 = dm[!duplicated(dm$id), c("id", "class", zpreds, "Gender")]; nd0$n = 1
w = as.numeric(table(dm$id)[nd0$id]); w = w / sum(w)  # weight children by number of valid items
for (m in intersect(c("cat_class", "cat_item"), names(fits))) {
  fit = fits[[m]]; set.seed(1); ids = sample(ndraws(fit), min(500, ndraws(fit)))
  ep = function(nd) { p = posterior_epred(fit, newdata = nd, re_formula = ~ (1 | id) + (1 | class), draw_ids = ids)
                      apply(p, c(1, 3), function(x) sum(x * w)) }  # draws x categories
  out = list()
  for (v in c(zpreds, "Gender")) {
    lo = nd0; hi = nd0
    if (v == "Gender") { lo$Gender = "F"; hi$Gender = "M" } else { lo[[v]] = lo[[v]] - 0.5; hi[[v]] = hi[[v]] + 0.5 }
    dd = ep(hi) - ep(lo)
    for (k in 1:3) out[[length(out) + 1]] = data.frame(model = m, term = v, category = dimnames(dd)[[2]][k], t(summ(dd[, k])))
  }
  ame = do.call(rbind, out); ame[, 4:7] = round(ame[, 4:7], 4); print(ame, row.names = FALSE)
  write.csv(ame, paste0("results/03_ame_", m, ".csv"), row.names = FALSE)
}
sink()
