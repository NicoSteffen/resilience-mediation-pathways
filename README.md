# Resilience as a mediator between life events and depressive symptoms

Analysis code for:

> Steffen, N., Sarubin, N., Reinhard, M. A., Jobst-Heel, A., Barton, B. B.,
> Kunz, J. I., Wolf, J., Bühner, M., Übleis, A., Padberg, F., & Goerigk, S.
> *Associations between adaptive life events and depressive symptoms are
> differentially mediated by resilience across the lifespan: path analyses in
> a naturalistic inpatient sample with major depression.*

## Contents

| File | Purpose |
|---|---|
| `00_functions.R` | Shared setup, constants, TAQ/CD-RISC scoring and table helpers |
| `01_main_discovery_analysis.R` | Discovery sample (N = 133 MDD inpatients) |
| `02_external_validation_analysis.R` | Validation sample (N = 67, PDD/BPD) |
| `03_joint_developmental_model.R` | Joint model over all three periods, with contrasts |
| `04_covariate_sensitivity_discovery.R` | Covariate-adjusted model, discovery sample |
| `05_covariate_sensitivity_validation.R` | Covariate-adjusted model, validation sample |
| `run_all.R` | Runs everything in the required order |

## Run order

All scripts are run from the project root and begin with
`source("00_functions.R")`.

```
01_main_discovery_analysis.R          independent
02_external_validation_analysis.R     independent
03_joint_developmental_model.R        requires 01
04_covariate_sensitivity_discovery.R  requires 01
05_covariate_sensitivity_validation.R requires 02
```

Scripts 03–05 reuse the data frames and fitted models that 01 and 02 leave in
the environment (`swdata`, `fit_lifetime`; `val`, `fit_validation`), so they
must run in the same session as the script they depend on. Each of them checks
for the objects it needs and stops with a message naming the required script,
rather than failing further down.

The simplest way to reproduce everything is:

```r
source("run_all.R")     # or: Rscript run_all.R
```

Runtime is dominated by the bootstrap (`N_BOOTSTRAP <- 10000`) and is roughly
half an hour on a current desktop machine.

Each script writes to its own folder:

| Script | Output folder |
|---|---|
| `01` | `outputs_main_discovery_analysis/` |
| `02` | `outputs_external_validation/` |
| `03` | `outputs_joint_developmental/` |
| `04` | `outputs_covariate_sensitivity/` |
| `05` | `outputs_covariate_sensitivity_validation/` |

## Data

Raw data are individual-level patient data and are **not published**. They are
available from the corresponding author on reasonable request, subject to the
approvals of the ethics committee of the Faculty of Medicine, LMU Munich
(Nr. 531-13).

Each analysis script names the SPSS export it reads in a `DATA_FILE` constant
at the top; both are expected under `raw_data/` relative to the project root.

## Software

R 4.5.0. Attached packages: `tidyverse`, `lavaan`. Called via `::`:
`foreign`, `psych`, `broom`, `Hmisc`, `report`, `flextable`.

The bootstrap seed is fixed in `00_functions.R` (`SEED <- 123`) and applied to
every model.

## Scoring

* **TAQ** (Hofmann et al., 1999), scored per the manual
  (Spinazzola et al., 2019). The two "family secrets" items and the check-up
  item 43 are excluded; the remaining 40 items are renumbered 1–40. The German
  version covers three developmental periods — childhood (7–12), youth (13–18),
  adulthood (>18) — whereas the English version additionally covers ages 0–6.
  Subscale scores are the mean of *endorsed* items only (a response of 0 is
  treated as "not endorsed"); a respondent endorsing nothing scores 0.
  TAQ⁺ / TAQ⁻ lifetime composites average the across-period subscale sums
  (range 0–9); period-specific composites average the period subscale means
  (range 0–3).
* **CD-RISC**: cases with more than `MAX_MISSING_CDRISC` (= 2) missing items
  are excluded; remaining missing items are person-mean imputed.
* **Mediation models**: `lavaan::sem()`, `std.ov = TRUE`, `fixed.x = FALSE`,
  FIML, 10,000 bootstrap replications.

---

## Reverse-scoring of TAQ item 2

TAQ item 2 is positively worded ("someone you could rely on") and is the only
item that is reverse-scored.

The **originally submitted** analysis selected the reverse-scored columns with
the pattern `"TAQ_2"`. That pattern matches `TAQ_2a`–`TAQ_2c` but **also**
`TAQ_20a` through `TAQ_29c`. Because the reverse-coding step ran first and
converted those columns to numeric, the subsequent standard-coding step could
no longer match them, so items 20–29 kept the reversed scoring.

Items 20–29 contribute to Witnessing (20, 22, 23, 24), Neglect (21, 27),
Alcohol/Drugs (25, 26) and Physical Abuse (28, 29) — four of the eight
subscales making up TAQ⁻. The error therefore affected TAQ⁻ (lifetime and all
three periods), its correlations, its α, and every path involving TAQ⁻. TAQ⁺
was unaffected.

The code in this repository reverse-scores item 2 only. The pattern is
anchored, `"^TAQ_2[abc]$"`, so it cannot pick up items 20–29 again. The
reverse-scored columns are printed at run time as a check.

## Differences from the originally submitted analysis

These were identified while preparing the code for release and are corrected
here. Numbers in the current output may therefore differ from the submitted
tables.

1. **TAQ reverse-scoring** — see the section above. Affects TAQ⁻ throughout.
2. **No seed in the submitted discovery run.** `01` did not set a seed before
   the bootstrap (`02` did). A seed is now set for every model via `SEED`;
   expect third-decimal differences from the submitted Table 3 and Table 4.
3. **Two grouping variables in the validation script.** The original filtered
   on `group` in one place and on `group_1` in another. It is now set once via
   `GROUP_VAR`, and the script asserts `n == 67`.
4. **Incomplete developmental periods.** Period totals were summed with
   `na.rm = TRUE`, which counted a skipped developmental period as an absence
   of events and produced artificially low lifetime scores for three
   respondents. Such cases are now missing at the composite level and handled
   by FIML.
5. **CD-RISC exclusion threshold.** The Methods state that participants with
   >10% missing CD-RISC values were excluded, which for the 25-item version is
   `MAX_MISSING_CDRISC <- 2`. The original code used 3. In these data both
   rules select the same cases (no respondent had exactly 3 missing items), so
   the sample is unchanged. The same constant is applied to the 10-item
   CD-RISC of the validation sample, where it is likewise inconsequential (no
   respondent had 2 missing items; the seven excluded cases had all 10
   missing).
6. **`BDI_score` vs. `BDI_Score`.** A relabelling step used `"BDI_score"`,
   which never matched the actual column name, so Table S5 carried the raw
   variable name.
7. **Table row selection.** The lifetime and validation tables removed rows by
   hard-coded index while the three period tables did not, so the tables did
   not have the same rows. `tidy_mediation()` now applies one rule everywhere:
   regression paths, the covariance between the two exogenous predictors, and
   the defined indirect and total effects. **Estimates are unchanged** — only
   which rows are printed.
8. **Correlation table row labels.** The adverse-scale rows were labelled
   `TAQ_PC`, `TAQ_PY`, `TAQ_PA`, duplicating the adaptive-scale labels. They
   are now `TAQ_NC`, `TAQ_NY`, `TAQ_NA`.
9. **Benjamini–Hochberg correction** was applied to the full symmetric
   correlation matrix, counting each pair twice. It is now applied to the
   unique pairs only. BH is invariant to exact duplication, so the adjusted
   p-values are identical.
10. **Helper functions.** `roundallnumerics()` used `df[, j]` and silently did
    nothing when passed a tibble, so two tables were never rounded;
    `p_labeller()` left p = .01 unstarred and coerced its vector to character
    mid-loop, so later comparisons were string comparisons. Both are replaced
    (`round_numeric()`, `p_label()`).
11. **Silent recoding failures.** Unmatched response labels previously became
    `NA` without notice; `recode_items()` now warns and lists them.
12. **Undeclared TAQ response codes in the validation sample.** Nine yes/no TAQ
    columns of that export carry the bare numeric codes `5` and `3` instead of
    a value label. Taken at face value they exceed the 0–3 range of the scored
    items. Per the value-label table of the export they follow the Likert
    scheme of the remaining TAQ items (1 = never, 2 = rarely,
    3 = occasionally, 4 = often, 5 = don't know), and they are now scored
    accordingly (`5` → 0, `3` → 2; reverse-scored `3` → 0).
