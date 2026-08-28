# ==============================================================================
# 05_covariate_sensitivity_validation.R
#
# Covariate-adjusted sensitivity analysis for the mediation model of the
# validation sample (Table S5), reported alongside the unadjusted estimates.
#
# Covariates: sex, age, diagnosis, medication. Age at onset is NOT included --
# it was not consistently recorded in this sample, which is why it is not
# reported in Table 1. Medication enters as the overall psychotropic indicator
# rather than the individual drug classes: with N = 67 the class-specific cells
# are thin.
#
# Produces, in outputs_covariate_sensitivity_validation/:
#   covariate_sensitivity_validation.docx / .csv
#   covariate_sensitivity_validation_log.txt
#
# Requires 02_external_validation_analysis.R to have run first
# (val, fit_validation, GROUP_VAR).
# Run from the project root. See README.md for the run order.
# ==============================================================================

source("00_functions.R")

require_objects("val", "fit_validation", "GROUP_VAR",
                after = "02_external_validation_analysis.R")

OUTPUT_DIR <- "outputs_covariate_sensitivity_validation"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Mediator and outcome
# ==============================================================================

# Set locally rather than inherited, and checked against the fitted main model.
# This sample uses the CD-RISC-10; the discovery sample uses the CD-RISC-25.
SENS_VAL_MEDIATOR <- "CDRISC_10_sum"
SENS_VAL_OUTCOME  <- "BDI_Score"

stopifnot(all(c(SENS_VAL_MEDIATOR, SENS_VAL_OUTCOME) %in% names(val)))
check_model_vars(fit_validation, SENS_VAL_MEDIATOR, SENS_VAL_OUTCOME,
                 what = "The covariate-adjusted validation model")


# ==============================================================================
# 2. Covariates
# ==============================================================================

# Sex: coded by value, not by factor order. The factor levels in this sample
# start with "männlich", the reverse of the discovery sample, so
# as.numeric(factor) - 1 would flip the coding. Table S1 codes 0 = female,
# 1 = male; this reproduces it in both samples.
val$Sex_num <- as.numeric(val$Geschlecht == "männlich")

# Diagnosis: 0 = CD (persistent depressive disorder), 1 = BPD.
val$Dx_bpd <- as.numeric(val[[GROUP_VAR]] == "BPD")

# Medication: any psychotropic prescribed at admission, 0 = no, 1 = yes.
val$Med_num <- as.numeric(val$medication == "yes")

COVARIATES_VAL <- c("Sex_num", "Alter_echt", "Dx_bpd", "Med_num")

# Verify the coding before trusting the model.
print(table(val$Geschlecht,      val$Sex_num, useNA = "ifany"))
print(table(val[[GROUP_VAR]],    val$Dx_bpd,  useNA = "ifany"))
print(table(val$medication,      val$Med_num, useNA = "ifany"))
print(colSums(is.na(val[, COVARIATES_VAL])))

# A covariate with no variance, or with a near-empty cell, will not support a
# bootstrap; stop early rather than produce unstable estimates.
stopifnot(all(vapply(val[COVARIATES_VAL], function(x) {
  var(x, na.rm = TRUE) > 0
}, logical(1))))


# ==============================================================================
# 3. Adjusted model
# ==============================================================================

# Same estimation settings as fit_mediation(). Covariates enter BOTH equations:
# adjusting only the outcome would leave the a path, and therefore the indirect
# effect, unadjusted. Exogenous variables covary freely, which keeps the model
# saturated (df = 0) as in the unadjusted model.

model_adjusted_val <- sprintf(
  "
  %1$s ~ a*TAQ_Pos + b*TAQ_Neg + Sex_num + Alter_echt + Dx_bpd + Med_num
  %2$s ~ c*%1$s
  %2$s ~ d*TAQ_Pos + e*TAQ_Neg + Sex_num + Alter_echt + Dx_bpd + Med_num

  TAQ_Pos    ~~ TAQ_Neg + Sex_num + Alter_echt + Dx_bpd + Med_num
  TAQ_Neg    ~~ Sex_num + Alter_echt + Dx_bpd + Med_num
  Sex_num    ~~ Alter_echt + Dx_bpd + Med_num
  Alter_echt ~~ Dx_bpd + Med_num
  Dx_bpd     ~~ Med_num

  ac        := a*c
  bc        := b*c
  total.pos := d + (a*c)
  total.neg := e + (b*c)
  ",
  SENS_VAL_MEDIATOR, SENS_VAL_OUTCOME
)

cat(model_adjusted_val)

set.seed(SEED)
fit_validation_adj <- lavaan::sem(
  model_adjusted_val,
  data      = val,
  fixed.x   = FALSE,
  std.ov    = TRUE,
  missing   = "fiml",
  se        = "bootstrap",
  bootstrap = N_BOOTSTRAP
)

print(summary(fit_validation_adj))

message("df unadjusted: ", lavaan::fitMeasures(fit_validation, "df"),
        " | df adjusted: ", lavaan::fitMeasures(fit_validation_adj, "df"))


# ==============================================================================
# 4. Side-by-side comparison
# ==============================================================================

# Longer names first: tidy_mediation() replaces labels sequentially.
labels_val <- c(stats::setNames(c("CDRISC", "BDI"),
                                c(SENS_VAL_MEDIATOR, SENS_VAL_OUTCOME)),
                TAQ_Pos = "TAQ+", TAQ_Neg = "TAQ-")

val_unadj <- tidy_mediation(fit_validation,     labels = labels_val)
val_adj   <- tidy_mediation(fit_validation_adj, labels = labels_val)

keep_cols <- c("Term", "Point estimate", "p-value", "CI lower", "CI upper")
val_unadj <- val_unadj[, keep_cols]
val_adj   <- val_adj[,   keep_cols]

names(val_unadj)[-1] <- paste0(names(val_unadj)[-1], " (unadjusted)")
names(val_adj)[-1]   <- paste0(names(val_adj)[-1],   " (adjusted)")

# all.x = TRUE keeps the rows of the unadjusted model and drops the covariate
# paths, which are not the point of the table.
sens_val <- merge(val_unadj, val_adj, by = "Term", all.x = TRUE, sort = FALSE)
sens_val <- sens_val[match(val_unadj$Term, sens_val$Term), ]
rownames(sens_val) <- NULL

print(sens_val)
save_table(sens_val, "covariate_sensitivity_validation.docx")
write.csv(sens_val,
          file.path(OUTPUT_DIR, "covariate_sensitivity_validation.csv"),
          row.names = FALSE)


# ==============================================================================
# 5. Log
# ==============================================================================

sens_log_val <- c(
  "Covariate-adjusted sensitivity analysis -- validation sample",
  "Covariates: sex, age, diagnosis (BPD vs CD), any psychotropic medication.",
  "Age at onset excluded: not consistently recorded in this sample.",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Seed:", SEED, "| Bootstrap replications:", N_BOOTSTRAP),
  "",
  "--- Covariate coding ---",
  capture.output(print(table(val$Geschlecht,   val$Sex_num, useNA = "ifany"))),
  capture.output(print(table(val[[GROUP_VAR]], val$Dx_bpd,  useNA = "ifany"))),
  capture.output(print(table(val$medication,   val$Med_num, useNA = "ifany"))),
  "",
  "--- Missing values on covariates ---",
  capture.output(print(colSums(is.na(val[, COVARIATES_VAL])))),
  "",
  "--- Degrees of freedom ---",
  paste("df unadjusted:", lavaan::fitMeasures(fit_validation, "df"),
        "| df adjusted:", lavaan::fitMeasures(fit_validation_adj, "df")),
  "",
  "--- Unadjusted vs. adjusted ---",
  capture.output(print(sens_val)),
  "",
  "--- Adjusted model, full summary ---",
  capture.output(summary(fit_validation_adj)),
  "",
  "--- Session info ---",
  capture.output(sessionInfo())
)

writeLines(sens_log_val,
           file.path(OUTPUT_DIR, "covariate_sensitivity_validation_log.txt"))
message("Written: ",
        file.path(OUTPUT_DIR, "covariate_sensitivity_validation_log.txt"))
