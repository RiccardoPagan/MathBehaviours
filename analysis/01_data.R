#### 01. DATA PREPARATION ####
# Input:  data/DatasetBehaviours.csv (long format: 32 rows per child at T0, 32 at T1)
# Output: analysis/results/trials.csv   one row per T1 trial (outcome + predictors)
#         analysis/results/children.csv one row per child (predictors, scale items)
# Run from the project root (open MathBehaviours.Rproj, or: Rscript analysis/01_data.R)

if (basename(getwd()) == "analysis") setwd("..")
d = read.csv("data/DatasetBehaviours.csv")
dir.create("analysis/results", showWarnings = FALSE)

# vo4a16 has its 32 T1 rows exported twice (exact copies): drop exact duplicate rows
cat("Exact duplicate rows removed:", sum(duplicated(d)), "\n")
d = d[!duplicated(d), ]

# choice: exactly one of skip/help/self is 1; otherwise missing
d$choice = NA
d$choice[d$skip %in% 1] = "skip"
d$choice[d$help %in% 1] = "help"
d$choice[d$self %in% 1] = "self"

# class = school + grade + section, from the code (e.g. "dv4a01" -> "dv4a"); the Class column is empty
d$class = substr(d$id, 1, 4)
d$school = substr(d$id, 1, 2)

#### T1 trials ####
d1 = d[d$session == 2, ]
# ASSUMPTION: rows are in presentation order within child (check with Riccardo)
d1$pos = ave(seq_along(d1$id), d1$id, FUN = seq_along)
d1$block = ifelse(d1$pos <= 16, 1, 2)

#### one row per child ####
preds = c("worry", "sc_state", "math_anxiety", "sc_trait", "buoyancy", "avoidance", "math_ability_T0")
items = c(paste0("trait_buoy_", 1:4), paste0("trait_avoidance_", 1:10))
ch = d1[!duplicated(d1$id), c("id", "class", "school", "Gender", "Grade", preds, items)]

# z-scores computed on children (not on trials), as in the preregistration
for (v in preds) ch[[paste0(v, "_z")]] = as.numeric(scale(ch[[v]]))

# per-child choice counts and proportions (missing choices excluded from the denominator)
ch$n_valid = as.numeric(tapply(!is.na(d1$choice), d1$id, sum)[ch$id])
ch$n_missing = as.numeric(tapply(is.na(d1$choice), d1$id, sum)[ch$id])
for (k in c("skip", "help", "self")) {
  ch[[paste0("n_", k)]] = as.numeric(tapply(d1$choice %in% k, d1$id, sum)[ch$id])
  ch[[paste0("p_", k)]] = ch[[paste0("n_", k)]] / ch$n_valid
}

#### merge predictors into trials ####
tr = merge(d1[, c("id", "item", "pos", "block", "accuracy", "choice")],
           ch[, c("id", "class", "school", "Gender", "Grade", preds, paste0(preds, "_z"))], by = "id")
tr = tr[order(tr$id, tr$pos), ]

#### checks ####
cat("Children:", nrow(ch), " classes:", length(unique(ch$class)), " schools:", length(unique(ch$school)), "\n")
cat("T1 rows per child:\n"); print(table(table(d1$id)))
cat("Children with T0 data:", length(unique(d$id[d$session == 1])), "\n")
cat("Missing values in predictors:", sum(is.na(ch[, preds])), "\n")
cat("T1 trials:", nrow(tr), " with valid choice:", sum(!is.na(tr$choice)), "\n")
print(table(tr$choice, useNA = "always"))

write.csv(tr, "analysis/results/trials.csv", row.names = FALSE)
write.csv(ch, "analysis/results/children.csv", row.names = FALSE)
