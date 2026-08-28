# ==============================================================================
# 01_main_discovery_analysis.R
#
# Discovery sample: N = 133 psychiatric inpatients with current MDD.
#
# Tests whether resilience (CD-RISC-25) mediates the association between
# adaptive (TAQ+) and adverse (TAQ-) life events and depressive symptoms
# (BDI-II), for lifetime scores and separately for childhood, youth and
# adulthood.
#
# Produces, in outputs_main_discovery_analysis/:
#   Table1_sample_discovery.docx     sample characteristics
#   Table1_sample_by_sex.docx        sample characteristics by sex
#   Reliability_discovery.docx       internal consistency (Table S4)
#   Table2_correlations.docx         descriptives and zero-order correlations
#   Table2_correlations_p_BH.docx    BH-adjusted p-values for Table 2
#   TableS1_covariates.docx          BDI-II on candidate covariates
#   TableS2_period_regressions.docx  CD-RISC / BDI-II on period-specific TAQ
#   Table3_mediation_lifetime.docx   lifetime mediation model
#   Table4a/b/c_mediation_*.docx     mediation models by developmental period
#
# First script in the pipeline; requires nothing but 00_functions.R.
# Run from the project root. See README.md for the run order.
# ==============================================================================

source("00_functions.R")

OUTPUT_DIR <- "outputs_main_discovery_analysis"

# SPSS export of the discovery sample. Not published (patient data);
# see README.md for access conditions.
DATA_FILE <- file.path("raw_data", "file.sav")


# ==============================================================================
# 1. Import
# ==============================================================================

swdata <- foreign::read.spss(DATA_FILE, to.data.frame = TRUE)


# ==============================================================================
# 2. Preprocessing
# ==============================================================================

## 2.1 CD-RISC completeness ----------------------------------------------------
# Cases with more than MAX_MISSING_CDRISC missing items are excluded; the
# remainder are person-mean imputed in section 4.

cdrisc_cols <- grep("^CDRISC_T0_", names(swdata), value = TRUE)
missing_cdrisc <- rowSums(is.na(swdata[, cdrisc_cols]))

message("CD-RISC missing items per case:")
print(table(missing_cdrisc))

swdata <- swdata[missing_cdrisc <= MAX_MISSING_CDRISC, ]


## 2.2 Variable names and factor levels ----------------------------------------

swdata <- dplyr::rename(
  swdata,
  ID               = Code,
  Gender           = Geschlecht,
  Age              = Alter,
  Weight           = Gewicht_kg,
  Height           = `Größe_cm`,
  Comorbidities    = Anzahl_Diagnosen,
  Age_onset        = Alter_bei_Ersterkrankung,
  Hospitalizations = Anzahl_Klinikaufenthalte,
  Medication       = Aufnahne_Medikation_ja_nein,
  Antidepressants  = Aufnahme_Antidepressiva,
  Antipsychotics   = Aufnahme_Antipsychotika,
  Benzodiazepines  = Aufnahme_Benzodiazepine
)

# Levels are relabelled by position, so the order below must match the order
# of the SPSS value labels.
levels(swdata$Gender) <- c("Female", "Male")

for (v in c("Medication", "Antidepressants", "Antipsychotics", "Benzodiazepines")) {
  levels(swdata[[v]]) <- c("yes", "no")
  swdata[[v]] <- relevel(swdata[[v]], ref = "no")
}


## 2.3 Drop administrative and unused columns ----------------------------------
# TAQ_T0_43 is a check-up item and is not part of any subscale.

swdata <- swdata[, !names(swdata) %in%
                   c("filter_.", "Uhrzeit", "Bemerkungen", "Datum", "TAQ_T0_43")]

# Instruments and follow-up time points not used in this paper.
drop_patterns <- c("MSWS", "RS25", "SES", "SB",
                   "BDI_T2", "BDI_T4", "BDI_T6", "BDI_T8")
dropped <- unique(unlist(lapply(drop_patterns, grep, x = names(swdata), value = TRUE)))
message("Dropping ", length(dropped), " unused columns:")
print(dropped)
swdata <- swdata[, !names(swdata) %in% dropped]

names(swdata) <- gsub("TAQ_T0", "TAQ", names(swdata))


## 2.4 Diagnostic groups and selection of the MDD subsample --------------------

swdata$Diagnosis <- dplyr::case_when(
  swdata$Hauptdiagnose_Detail == "Bipolar" ~ "Bipolar",
  swdata$Hauptdiagnose_ICD10_Kategorie ==
    "Schizophrenie, schizotype und wahnhafte Störungen (F20-F29)" ~ "Psychotic",
  swdata$Hauptdiagnose_ICD10_Kategorie ==
    "Affektive Störungen (F30-F39)" ~ "Depression",
  swdata$Hauptdiagnose_ICD10_Kategorie ==
    "Persönlichkeits- und Verhaltensstörungen (F60-F69)" ~ "PD",
  swdata$Hauptdiagnose_ICD10_Kategorie ==
    "Neurotische, Belastungs- und somatoforme Störungen (F40-F49)" ~
    "Neurotic, stress-related and somatoform",
  TRUE ~ NA_character_
)

# Sample characteristics reported throughout
characteristics <- c("Gender", "Age", "Comorbidities", "Age_onset",
                     "Hospitalizations", "Medication", "Antidepressants",
                     "Antipsychotics", "Benzodiazepines")

message("Full cohort by diagnostic group:")
print(as.data.frame(
  report::report_sample(swdata[, c(characteristics, "Diagnosis")],
                        group_by = "Diagnosis")
))

# Analysis sample: unipolar MDD only (bipolar I and dysthymia excluded).
swdata <- swdata[
  swdata$Hauptdiagnose_ICD10_Kategorie == "Affektive Störungen (F30-F39)" &
    !swdata$Hauptdiagnose_Detail %in% c("Bipolar", "Dystymie"),
]

message("Analysis sample: n = ", nrow(swdata))
stopifnot(nrow(swdata) == 133)


## 2.5 TAQ item selection and renumbering --------------------------------------
# Items 9 and 24 ("family secrets") were dropped from later TAQ versions.
# The remaining 40 items x 3 developmental periods are renumbered 1a-40c.

swdata <- swdata[, !grepl("TAQ_9|TAQ_24", names(swdata))]

taq_idx <- grep("TAQ", names(swdata))
stopifnot(length(taq_idx) == 120)   # 40 items x 3 periods
names(swdata)[taq_idx] <- paste0("TAQ_", rep(1:40, each = 3),
                                 rep(letters[1:3], times = 40))


# ==============================================================================
# 3. TAQ scoring
# ==============================================================================

TAQ_MAP <- c(
  "Oft?/sehr stark"               = 3,
  "Gelegentlich?/im Grunde schon" = 2,
  "Selten?/ein wenig"             = 0,
  "Nie?/gar nicht"                = 0,
  "ja"                            = 3,
  "nein"                          = 0,
  "Weiß nicht"                    = 0
)

TAQ_MAP_REVERSED <- c(
  "Oft?/sehr stark"               = 0,
  "Gelegentlich?/im Grunde schon" = 0,
  "Selten?/ein wenig"             = 2,
  "Nie?/gar nicht"                = 3,
  "Weiß nicht"                    = 0
)

# TAQ item 2 is positively worded and is the only reverse-scored item. The
# anchored pattern matters: "^TAQ_2" would also match items 20-29 and reverse
# four of the eight TAQ- subscales. See README.md, "Reverse-scoring of TAQ
# item 2".
reverse_cols <- grep("^TAQ_2[abc]$", names(swdata), value = TRUE)
message("Reverse-scored TAQ columns (n = ", length(reverse_cols), "): ",
        paste(reverse_cols, collapse = ", "))

standard_cols <- setdiff(grep("^TAQ_[0-9]", names(swdata), value = TRUE),
                         reverse_cols)

# Reverse-scored columns first: recoding turns labels into numbers, so a column
# recoded once can no longer be matched by the second map.
swdata <- recode_items(swdata, reverse_cols,  TAQ_MAP_REVERSED)
swdata <- recode_items(swdata, standard_cols, TAQ_MAP)

# Subscales, period totals and TAQ+ / TAQ- composites
swdata <- add_taq_scores(swdata)


## 3.1 Completeness of the TAQ composites --------------------------------------
# Reported in the Methods and in the notes to Tables 2, 3 and S2.

taq_n <- data.frame(
  Score = c("TAQ_Pos / TAQ_Neg (lifetime)",
            "TAQ_Pos_a / TAQ_Neg_a (childhood)",
            "TAQ_Pos_b / TAQ_Neg_b (youth)",
            "TAQ_Pos_c / TAQ_Neg_c (adulthood)"),
  n = c(sum(!is.na(swdata$TAQ_Pos)),
        sum(!is.na(swdata$TAQ_Pos_a)),
        sum(!is.na(swdata$TAQ_Pos_b)),
        sum(!is.na(swdata$TAQ_Pos_c)))
)
print(taq_n)


# ==============================================================================
# 4. CD-RISC and BDI-II scoring
# ==============================================================================

CDRISC_MAP <- c(
  "überhaupt nicht wahr" = 0,
  "selten wahr"          = 1,
  "manchmal wahr"        = 2,
  "oft wahr"             = 3,
  "fast immer wahr"      = 4
)

cdrisc_cols <- grep("^CDRISC_T0_", names(swdata), value = TRUE)
swdata <- recode_items(swdata, cdrisc_cols, CDRISC_MAP)
swdata <- impute_person_mean(swdata, cdrisc_cols)

swdata$CDRISC_25_T0_sum <- rowSums(swdata[, cdrisc_cols])

# CD-RISC-10 subset, for comparability with the validation sample. Reported in
# Table 1 only; the models below use the CD-RISC-25.
cdrisc10_items <- paste0("CDRISC_T0_", c(1, 4, 6, 7, 8, 11, 14, 16, 17, 19))
swdata$CDRISC_10_sum <- rowSums(swdata[, cdrisc10_items])


# ==============================================================================
# 5. Table 1 -- sample characteristics
# ==============================================================================

n_total <- nrow(swdata)

table1_discovery <- tibble::tibble(
  Variable = c(
    paste("Sex:",             names(n_percent(swdata$Gender))),
    "Major depression, recurrent",
    "Major depression, single episode",
    "Age (years)",
    "Age at onset (years)",
    "Comorbidities",
    "Hospitalizations",
    paste("Medication:",      names(n_percent(swdata$Medication))),
    paste("Antidepressants:", names(n_percent(swdata$Antidepressants))),
    paste("Antipsychotics:",  names(n_percent(swdata$Antipsychotics))),
    paste("Benzodiazepines:", names(n_percent(swdata$Benzodiazepines))),
    "TAQ+ (adaptive)",
    "TAQ- (adverse)",
    "CD-RISC-25",
    "CD-RISC-10",
    "BDI-II"
  ),
  Value = c(
    n_percent(swdata$Gender),
    sprintf("%d (%.1f%%)",
            sum(swdata$Hauptdiagnose_Detail == "Rezidivierend Depressiv", na.rm = TRUE),
            100 * mean(swdata$Hauptdiagnose_Detail == "Rezidivierend Depressiv", na.rm = TRUE)),
    sprintf("%d (%.1f%%)",
            sum(swdata$Hauptdiagnose_Detail == "Episodisch Depressiv (MD)", na.rm = TRUE),
            100 * mean(swdata$Hauptdiagnose_Detail == "Episodisch Depressiv (MD)", na.rm = TRUE)),
    mean_sd(swdata$Age),
    mean_sd(swdata$Age_onset),
    mean_sd(swdata$Comorbidities),
    mean_sd(swdata$Hospitalizations),
    n_percent(swdata$Medication),
    n_percent(swdata$Antidepressants),
    n_percent(swdata$Antipsychotics),
    n_percent(swdata$Benzodiazepines),
    mean_sd(swdata$TAQ_Pos),
    mean_sd(swdata$TAQ_Neg),
    mean_sd(swdata$CDRISC_25_T0_sum),
    mean_sd(swdata$CDRISC_10_sum),
    mean_sd(swdata$BDI_T0_Gesamt)
  )
)
names(table1_discovery)[2] <- paste0("Discovery sample (n = ", n_total, ")")

print(table1_discovery, n = Inf)
save_table(table1_discovery, "Table1_sample_discovery.docx")

# Breakdown by sex, for the table footnote
table1_by_sex <- as.data.frame(
  report::report_sample(swdata[, c(characteristics)], by = "Gender")
)
save_table(table1_by_sex, "Table1_sample_by_sex.docx")


# ==============================================================================
# 6. Internal consistency (Table S4)
# ==============================================================================

items_bdi <- grep("^BDI_T0", names(swdata), value = TRUE)
items_bdi <- setdiff(items_bdi, c("BDI_T0_S_a", "BDI_T0_S_a_Aus",
                                  "BDI_T0_S_Aus", "BDI_T0_Gesamt"))
swdata[items_bdi] <- lapply(swdata[items_bdi], \(x) as.numeric(as.character(x)))

reliability <- reliability_table(swdata, cdrisc_cols, items_bdi)
print(reliability, n = Inf)
save_table(reliability, "Reliability_discovery.docx")


# ==============================================================================
# 7. Table 2 -- descriptives and zero-order correlations
# ==============================================================================

cor_vars <- c("TAQ_Pos", "TAQ_Pos_a", "TAQ_Pos_b", "TAQ_Pos_c",
              "TAQ_Neg", "TAQ_Neg_a", "TAQ_Neg_b", "TAQ_Neg_c",
              "CDRISC_25_T0_sum", "BDI_T0_Gesamt")

cor_labels <- c("TAQ_P", "TAQ_PC", "TAQ_PY", "TAQ_PA",
                "TAQ_N", "TAQ_NC", "TAQ_NY", "TAQ_NA",
                "CDRISC", "BDI")

cor_data <- swdata[, cor_vars]

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
dimnames(cor_matrix) <- list(cor_labels, cor_labels)

table2 <- data.frame(
  Variable = cor_labels,
  M  = sapply(cor_data, mean, na.rm = TRUE),
  SD = sapply(cor_data, sd,   na.rm = TRUE),
  cor_matrix,
  row.names = NULL,
  check.names = FALSE
)
table2 <- round_numeric(table2, 2)
save_table(table2, "Table2_correlations.docx")

# Benjamini-Hochberg FDR correction over the unique variable pairs. Applying it
# to the full symmetric matrix would count every pair twice.
p_matrix <- Hmisc::rcorr(as.matrix(cor_data))$P
dimnames(p_matrix) <- list(cor_labels, cor_labels)

lower <- lower.tri(p_matrix)
p_adj <- p_matrix
p_adj[] <- NA_real_
p_adj[lower] <- p.adjust(p_matrix[lower], method = "BH")
p_adj[upper.tri(p_adj)] <- t(p_adj)[upper.tri(p_adj)]

table2_p <- data.frame(Variable = cor_labels, round(p_adj, 3),
                       row.names = NULL, check.names = FALSE)
save_table(table2_p, "Table2_correlations_p_BH.docx")


# ==============================================================================
# 8. Table S1 -- candidate covariates
# ==============================================================================

covariate_model <- lm(
  BDI_T0_Gesamt ~ Gender + Age + Age_onset + CDRISC_25_T0_sum,
  data = swdata
)
print(summary(covariate_model))

tableS1 <- round_numeric(broom::tidy(covariate_model, conf.int = TRUE), 3)
save_table(tableS1, "TableS1_covariates.docx")


# ==============================================================================
# 9. Table S2 -- period-specific regressions
# ==============================================================================

predictors <- c("TAQ_Pos_a", "TAQ_Pos_b", "TAQ_Pos_c",
                "TAQ_Neg_a", "TAQ_Neg_b", "TAQ_Neg_c")
outcomes   <- c("CDRISC_25_T0_sum", "BDI_T0_Gesamt")

simple_regression <- function(outcome, predictor, data) {
  fit <- lm(as.formula(paste(outcome, "~", predictor)), data = data)
  broom::tidy(fit, conf.int = TRUE) |>
    dplyr::filter(term == predictor) |>
    dplyr::transmute(
      Regression        = paste(outcome, "\u2190", predictor),
      `Point estimate`  = estimate,
      `Standard error`  = std.error,
      `t-value`         = statistic,
      `p-value`         = p.value,
      `CI lower`        = conf.low,
      `CI upper`        = conf.high
    )
}

tableS2 <- tidyr::expand_grid(outcome = outcomes, predictor = predictors) |>
  purrr::pmap(\(outcome, predictor) simple_regression(outcome, predictor, swdata)) |>
  dplyr::bind_rows() |>
  round_numeric(3)

print(tableS2)
save_table(tableS2, "TableS2_period_regressions.docx")


# ==============================================================================
# 10. Tables 3 and 4 -- mediation models
# ==============================================================================

# Mediator and outcome of this sample. Set here and passed as arguments; the
# follow-up scripts define their own and check them against fit_lifetime.
DISC_MEDIATOR <- "CDRISC_25_T0_sum"
DISC_OUTCOME  <- "BDI_T0_Gesamt"

stopifnot(all(c(DISC_MEDIATOR, DISC_OUTCOME) %in% names(swdata)))

base_labels <- c(CDRISC_25_T0_sum = "CDRISC", BDI_T0_Gesamt = "BDI")

## 10.1 Table 3 -- lifetime model ----------------------------------------------

fit_lifetime <- fit_mediation(swdata, "TAQ_Pos", "TAQ_Neg",
                              DISC_MEDIATOR, DISC_OUTCOME)
print(summary(fit_lifetime))

table3 <- tidy_mediation(
  fit_lifetime,
  labels = c(base_labels, TAQ_Pos = "TAQ+", TAQ_Neg = "TAQ-")
)
print(table3)
save_table(table3, "Table3_mediation_lifetime.docx")


## 10.2 Table 4 -- models by developmental period ------------------------------

periods <- list(
  a = list(label = "Childhood", file = "Table4a_mediation_childhood.docx"),
  b = list(label = "Youth",     file = "Table4b_mediation_youth.docx"),
  c = list(label = "Adulthood", file = "Table4c_mediation_adulthood.docx")
)

fits_period   <- list()
tables_period <- list()

for (p in names(periods)) {
  x_pos <- paste0("TAQ_Pos_", p)
  x_neg <- paste0("TAQ_Neg_", p)

  fits_period[[p]] <- fit_mediation(swdata, x_pos, x_neg,
                                    DISC_MEDIATOR, DISC_OUTCOME)
  message("\n--- ", periods[[p]]$label, " ---")
  print(summary(fits_period[[p]]))

  # Longer names first: tidy_mediation() replaces labels sequentially.
  labels <- c(base_labels,
              stats::setNames(
                c(sprintf("TAQ+ (%s)", tolower(periods[[p]]$label)),
                  sprintf("TAQ- (%s)", tolower(periods[[p]]$label))),
                c(x_pos, x_neg)
              ))

  tables_period[[p]] <- tidy_mediation(fits_period[[p]], labels = labels)
  save_table(tables_period[[p]], periods[[p]]$file)
}


# ==============================================================================
# 11. Session information
# ==============================================================================

sessionInfo()
