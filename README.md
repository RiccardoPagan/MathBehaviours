# Math Behaviours

How math attitudes (trait and state) relate to children's item-by-item choices in a math task: skip the problem, ask for a hint, or solve it alone. Grades 4 and 5, two sessions (T0, T1), N = 547.

## Folders

| Folder | Content |
|---|---|
| `data/` | Real data: `DatasetBehaviours.csv` (long format, 32 trials per child per session) |
| `analysis/` | **Final analysis scripts** (`01`–`05`), their outputs in `analysis/results/` and `analysis/figures/` |
| `Report/` | Manuscript (`paper.qmd`, apaquarto) and working notes (`note.md`) |
| `Slide/` | AIRIPA 2026 talk |
| `Preregistration/` | Frozen preregistration material: OSF text, earlier draft, simulated data, simulation and power-analysis scripts, planned analysis scripts, project notes (Nov 2025) |
| `archive/` | Superseded material, kept for reference only: first run of the analyses on the real data (`first_real_data_run/`), preliminary analyses of 31/10/2026 (`prelim_analyses_2026-10-31/`), a script from another project (`CheckDatiFeedback.R`) |

## How to run the final analyses

Open `MathBehaviours.Rproj` and run the scripts **in order**, from the project root or `analysis/`. Each script prints its results to the console and also saves them in `analysis/results/` (a `.txt` with the full printout, plus compact tables and numerical summaries).

Packages: `install.packages(c("brms", "lme4", "lavaan", "semTools", "psych", "loo", "posterior", "ggplot2"))`. Tested with R 4.5.1, brms 2.22.0.

| Script | What it does | Paper section | Time |
|---|---|---|---|
| `01_data.R` | Cleans the data (removes the duplicated export of `vo4a16`), codes the choice, computes z-scores on children, writes `trials.csv` and `children.csv` | — | seconds |
| `02_descriptives.R` | Sample, descriptives, correlations, choice distribution, accuracy, CFA of buoyancy and avoidance, bivariate table, VIF | Descriptive Statistics, Psychometric Screening, Assumption Check (Table of correlations with proportions) | ~1 min |
| `03_preregistered.R` | Preregistered ordinal pipeline: screening → M0/M1/M2 → category-specific effects → one-tailed *p* + FDR → indirect effects | Ordinal Models | 1–3 h |
| `04_multinomial.R` | Cumulative vs multinomial vs adjacent-category (LOO); multinomial model with all three contrasts; class random intercept; figure of predicted probabilities | Assumption Check, Multinomial Model, Robustness (class) | 1–3 h |
| `05_robustness.R` | Linear models on per-child scores; binary mixed models (self vs not, skip vs help, skip vs not), full and one attitude at a time | Robustness | ~5 min |

`03` and `04` use 4 cores each and are independent: they can run at the same time in two R sessions.

### Before running `03_preregistered.R`

Three decisions are set at the top of the script (section `DECISIONS TO SET BEFORE RUNNING`):

1. **Which CFI/RMSEA** the preregistered criterion (CFI > .95, RMSEA < .08) refers to. With the scaled (robust) indices avoidance fails (CFI = .92) and is excluded from the ordinal models; with the unscaled ones it passes (CFI = .97).
2. **Correlated residuals for avoidance** (at most two, only between items with clear overlap in content or wording). Candidates are listed in `02_descriptives.txt`; check the item wording before adding any.
3. **Order of the category-specific extensions** in Stage 2, which the preregistration does not fix (default: worry, perceived competence, buoyancy, avoidance).

### Notes

- Bayesian models are fitted on counts per child × choice (with `weights()`): the likelihood is identical to the trial-level one and sampling is about 15 times faster. LOO is then reconstructed at the trial level, as preregistered.
- `03` and `04` keep fitted models and MCMC draws **in memory only**. They save no fit caches or chain files in the project. Each run therefore fits the models again. `MathBehaviours.Rproj` also disables automatic workspace saving/restoration, so exiting RStudio does not save those fits in `.RData`.
- Both scripts automatically save global parameter summaries and diagnostics (`03_parameters.csv`, `03_diagnostics.csv`, `04_parameters.csv`, `04_diagnostics.csv`). The TXT printouts contain model specifications and priors. `03_spec.rds` records screening decisions, selected terms and sampling settings, without fitted models. The existing coefficient, contrast, LOO, indirect-effect and prediction tables remain the report inputs. `02_report.rds` contains only sample, choice, accuracy and psychometric summaries; `05` also saves model-fit statistics, random-effect SDs and convergence summaries. No separate compaction step is needed.
- Full MCMC caches were removed on 31/10/2026 after exporting and checking these summaries. All existing CSV/TXT outputs, LOO objects/comparisons and prediction summaries/figures were retained. New contrasts, predictions or sampling diagnostics beyond the saved summaries require refitting. Preregistration files, including the large `Results_vcovmtx.RData` workspace, were retained.
- `Report/paper.qmd` reads only compact final summaries from `analysis/results/` and the prediction figure from `analysis/figures/` (N = 547); all result tables are generated from these inputs. It uses no archived preliminary results or fitted models. The current ordinal run uses scaled CFA indices, no correlated residuals, and the extension order listed above; avoidance is excluded from that pipeline and retained only in exploratory follow-up models. Render from `Report/` with `quarto render paper.qmd --to all`.
- `01_data.R` assumes that rows are in presentation order within each child (used only for the block variable): to be confirmed.
- MCMC settings: 4 chains × 2,000 iterations (1,000 warmup), brms default priors, seed 2026.
