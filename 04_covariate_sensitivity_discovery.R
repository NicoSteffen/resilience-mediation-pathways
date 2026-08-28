# ==============================================================================
# 04_covariate_sensitivity_discovery.R
#
# Covariate-adjusted sensitivity analysis for the lifetime mediation model of
# the discovery sample (Table 3), reported alongside the unadjusted estimates.
#
# Covariates:
#   demographics     sex, age
#   illness history  age at onset, recurrent vs. single episode,
#                    number of comorbidities, number of hospitalisations
#   medication       antidepressants, antipsychotics, benzodiazepines
#
# Diagnosis is not a covariate here: the discovery sample is unipolar MDD
# throughout. The recurrent / single-episode distinction is its analogue and is
# included under illness history. Medication enters as the three drug classes
# rather than the overall yes/no indicator, which is close to constant in an
# inpatient sample and would carry little information.
#
# Produces, in outputs_covariate_sensitivity/:
#   TableS7_covariate_sensitivity.docx / .csv
#   covariate_sensitivity_log.txt
#
# Requires 01_main_discovery_analysis.R to have run first (swdata, fit_lifetime).
# Run from the project root. See README.md for the run order.
# ==============================================================================

source("00_functions.R")

require_objects("swdata", "fit_lifetime", after = "01_main_discovery_analysis.R")

OUTPUT_DIR <- "outputs_covariate_sensitivity"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Mediator and outcome
# ==============================================================================

# Set locally rather than inherited, and checked against the fitted main model:
# the adjusted model must use the same mediator and outcome as Table 3.
SENS_MEDIATOR <- "CDRISC_25_T0_sum"
SENS_OUTCOME  <- "BDI_T0_Gesamt"

stopifnot(all(c(SENS_MEDIATOR, SENS_OUTCOME) %in% names(swdata)))
check_model_vars(fit_lifetime, SENS_MEDIATOR, SENS_OUTCOME,
                 what = "The covariate-adjusted discovery model")


# ==============================================================================
# 2. Covariates
# ==============================================================================

# Coded by value rather than by factor order, so that the coding does not depend
# on how the SPSS labels happen to be ordered. Sex is 0 = female, 1 = male,
# matching Table S1.
swdata$Sex_num   <- as.numeric(swdata$Gender == "Male")
swdata$Recurrent <- as.numeric(
  swdata$Hauptdiagnose_Detail == "Rezidivierend Depressiv"
)
swdata$AD_num  <- as.numeric(swdata$Antidepressants  == "yes")
swdata$AP_num  <- as.numeric(swdata$Antipsychotics   == "yes")
swdata$BZD_num <- as.numeric(swdata$Benzodiazepines  == "yes")

swdata$Comorbid_n <- as.numeric(as.character(swdata$Comorbidities))

# Hospitalisations are right-skewed and are therefore log-transformed.
swdata$Hosp_log <- log1p(as.numeric(as.character(swdata$Hospitalizations)))

COVARIATES <- c("Sex_num", "Age", "Age_onset", "Recurrent",
                "Comorbid_n", "Hosp_log", "AD_num", "AP_num", "BZD_num")

# Verify the coding before trusting the model.
print(table(swdata$Gender,               swdata$Sex_num,   useNA = "ifany"))
print(table(swdata$Hauptdiagnose_Detail, swdata$Recurrent, useNA = "ifany"))
print(table(swdata$Antidepressants,      swdata$AD_num,    useNA = "ifany"))
print(table(swdata$Antipsychotics,       swdata$AP_num,    useNA = "ifany"))
print(table(swdata$Benzodiazepines,      swdata$BZD_num,   useNA = "ifany"))
print(summary(swdata[, c("Age", "Age_onset", "Comorbid_n", "Hosp_log")]))
print(colSums(is.na(swdata[, COVARIATES])))

# No covariate may be constant, or the model is not estimable.
stopifnot(all(vapply(swdata[COVARIATES],
                     \(x) var(x, na.rm = TRUE) > 0, logical(1))))

# Collinearity among the covariates. Age and age at onset correlate at about
# .70 in this sample; anything higher should be inspected before it is carried
# into the model.
cov_cor <- cor(swdata[, COVARIATES], use = "pairwise.complete.obs")
print(round(cov_cor, 2))


# ==============================================================================
# 3. Adjusted model
# ==============================================================================

# Same estimation settings as fit_mediation(). Covariates enter BOTH equations:
# adjusting only the outcome would leave the a path, and therefore the indirect
# effect, unadjusted. Exogenous variables covary freely, which keeps the model
# saturated (df = 0) as in the unadjusted model -- 55 pairs, so the covariance
# block is generated rather than typed out.

exog      <- c("TAQ_Pos", "TAQ_Neg", COVARIATES)
cov_pairs <- utils::combn(exog, 2)
cov_lines <- paste0("  ", cov_pairs[1, ], " ~~ ", cov_pairs[2, ],
                    collapse = "\n")

cov_terms <- paste(COVARIATES, collapse = " + ")

model_adjusted <- sprintf(
  "
  %1$s ~ a*TAQ_Pos + b*TAQ_Neg + %3$s
  %2$s ~ c*%1$s
  %2$s ~ d*TAQ_Pos + e*TAQ_Neg + %3$s

%4$s

  ac        := a*c
  bc        := b*c
  total.pos := d + (a*c)
  total.neg := e + (b*c)
  ",
  SENS_MEDIATOR, SENS_OUTCOME, cov_terms, cov_lines
)

cat(model_adjusted)

set.seed(SEED)
fit_lifetime_adj <- lavaan::sem(
  model_adjusted,
  data      = swdata,
  fixed.x   = FALSE,
  std.ov    = TRUE,
  missing   = "fiml",
  se        = "bootstrap",
  bootstrap = N_BOOTSTRAP
)

print(summary(fit_lifetime_adj))

message("df unadjusted: ", lavaan::fitMeasures(fit_lifetime, "df"),
        " | df adjusted: ", lavaan::fitMeasures(fit_lifetime_adj, "df"))

# Failed bootstrap draws are the first sign that the model is overloaded for
# this sample size.
message("Successful bootstrap draws: ",
        sum(!is.na(lavaan::lavInspect(fit_lifetime_adj, "boot")[, 1])),
        " of ", N_BOOTSTRAP)


# ==============================================================================
# 4. Table S7 -- side-by-side comparison
# ==============================================================================

# Longer names first: tidy_mediation() replaces labels sequentially.
labels_s7 <- c(stats::setNames(c("CDRISC", "BDI"), c(SENS_MEDIATOR, SENS_OUTCOME)),
               TAQ_Pos = "TAQ+", TAQ_Neg = "TAQ-")

s7_unadj <- tidy_mediation(fit_lifetime,     labels = labels_s7)
s7_adj   <- tidy_mediation(fit_lifetime_adj, labels = labels_s7)

keep_cols <- c("Term", "Point estimate", "p-value", "CI lower", "CI upper")
s7_unadj <- s7_unadj[, keep_cols]
s7_adj   <- s7_adj[,   keep_cols]

names(s7_unadj)[-1] <- paste0(names(s7_unadj)[-1], " (unadjusted)")
names(s7_adj)[-1]   <- paste0(names(s7_adj)[-1],   " (adjusted)")

# all.x = TRUE keeps the rows of the unadjusted model and drops the covariate
# paths, which are not the point of the table.
tableS7 <- merge(s7_unadj, s7_adj, by = "Term", all.x = TRUE, sort = FALSE)
tableS7 <- tableS7[match(s7_unadj$Term, tableS7$Term), ]
rownames(tableS7) <- NULL

print(tableS7)
save_table(tableS7, "TableS7_covariate_sensitivity.docx")
write.csv(tableS7, file.path(OUTPUT_DIR, "TableS7_covariate_sensitivity.csv"),
          row.names = FALSE)


# ==============================================================================
# 5. Log
# ==============================================================================

sens_log <- c(
  "Covariate-adjusted sensitivity analysis -- discovery sample",
  "Covariates: sex, age, age at onset, recurrent vs. single episode,",
  "comorbidities, log hospitalisations, antidepressants, antipsychotics,",
  "benzodiazepines.",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Seed:", SEED, "| Bootstrap replications:", N_BOOTSTRAP),
  "",
  "--- Covariate coding ---",
  capture.output(print(table(swdata$Gender, swdata$Sex_num, useNA = "ifany"))),
  capture.output(print(table(swdata$Hauptdiagnose_Detail, swdata$Recurrent,
                             useNA = "ifany"))),
  capture.output(print(table(swdata$Antidepressants, swdata$AD_num,
                             useNA = "ifany"))),
  capture.output(print(table(swdata$Antipsychotics, swdata$AP_num,
                             useNA = "ifany"))),
  capture.output(print(table(swdata$Benzodiazepines, swdata$BZD_num,
                             useNA = "ifany"))),
  "",
  "--- Continuous covariates ---",
  capture.output(print(summary(swdata[, c("Age", "Age_onset",
                                          "Comorbid_n", "Hosp_log")]))),
  "",
  "--- Missing values on covariates ---",
  capture.output(print(colSums(is.na(swdata[, COVARIATES])))),
  "",
  "--- Correlations among covariates ---",
  capture.output(print(round(cov_cor, 2))),
  "",
  "--- Degrees of freedom ---",
  paste("df unadjusted:", lavaan::fitMeasures(fit_lifetime, "df"),
        "| df adjusted:", lavaan::fitMeasures(fit_lifetime_adj, "df")),
  "",
  "--- Table S7: unadjusted vs. adjusted ---",
  capture.output(print(tableS7)),
  "",
  "--- Adjusted model, full summary ---",
  capture.output(summary(fit_lifetime_adj)),
  "",
  "--- Session info ---",
  capture.output(sessionInfo())
)

writeLines(sens_log, file.path(OUTPUT_DIR, "covariate_sensitivity_log.txt"))
message("Written: ", file.path(OUTPUT_DIR, "covariate_sensitivity_log.txt"))
