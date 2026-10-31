# Data preparation shared by all preliminary analyses
# (sourced by the other scripts; preliminary run, to be re-run later with more iterations)

d = read.csv("../../Analisys_truedata/DatasetBehaviours.csv")

d$choice = ifelse(d$skip %in% 1, "skip", ifelse(d$help %in% 1, "help", ifelse(d$self %in% 1, "self", NA)))
d$class = substr(d$id, 1, 4)  # e.g. "dv4a" = school + grade + section (Class column is all NA)

# vo4a16 has 64 rows at T1 (32 at T0): ambiguous, excluded for now pending clarification
d = subset(d, id != "vo4a16")

t0 = subset(d, session == 1)
d1 = subset(d, session == 2)
d1$pos = ave(seq_along(d1$id), d1$id, FUN = seq_along)  # presentation order assumed = row order
d1$block = ifelse(d1$pos <= 16, 1, 2)

preds = c("worry", "sc_state", "math_anxiety", "sc_trait", "buoyancy", "avoidance", "math_ability_T0")

# z-scores computed at the participant level (one row per child), then merged back
u = d1[!duplicated(d1$id), c("id", "class", "Gender", "Grade", preds)]
for (v in preds) u[[paste0(v, "_z")]] = as.numeric(scale(u[[v]]))
d1 = merge(d1[, c("id", "item", "pos", "block", "accuracy", "choice")], u, by = "id")
d1 = d1[order(d1$id, d1$pos), ]
zpreds = paste0(preds, "_z")

dm = subset(d1, !is.na(choice))  # analysis dataset (choice not missing)
dm$choice_c = factor(dm$choice, levels = c("skip", "help", "self"))
dm$choice_ord = ordered(dm$choice, levels = c("skip", "help", "self"))
dm$choice_seq = ordered(dm$choice, levels = c("self", "help", "skip"))  # sequential: self vs not, then help vs skip
dm$item = factor(dm$item)

rhs = paste(c(zpreds, "Gender"), collapse = " + ")
rhs_cs = paste(c(paste0("cs(", zpreds, ")"), "cs(Gender)"), collapse = " + ")
