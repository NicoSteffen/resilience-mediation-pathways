# ==============================================================================
# 03_joint_developmental_model.R
#
# Joint developmental model with contrasts between periods (discovery sample).
#
# The three period-specific models in Table 4 are estimated separately, so a
# significant effect in one period and a non-significant effect in another does
# not establish that the two differ. This script estimates ONE model containing
# all six period scores simultaneously and contrasts the period-specific
# indirect effects directly.
#
# Interpretation: in the joint model each path is the unique contribution of
# that period with the other five held constant. These estimates are therefore
# not comparable in size to those in Table 4, where each period was modelled on
# its own.
#
# Produces, in outputs_joint_developmental/:
#   TableS9_joint_developmental.docx / .csv
#   joint_developmental_log.txt
#
# Requires 01_main_discovery_analysis.R to have run first (swdata, fit_lifetime).
# Run from the project root. See README.md for the run order.
# ==============================================================================

source("00_functions.R")

require_objects("swdata", "fit_lifetime", after = "01_main_discovery_analysis.R")

OUTPUT_DIR <- "outputs_joint_developmental"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Mediator and outcome
# ==============================================================================

# Set locally rather than inherited, and checked against the fitted main model:
# the joint model must use the same mediator and outcome as Tables 3 and 4.
JOINT_MEDIATOR <- "CDRISC_25_T0_sum"
JOINT_OUTCOME  <- "BDI_T0_Gesamt"

stopifnot(all(c(JOINT_MEDIATOR, JOINT_OUTCOME) %in% names(swdata)))
check_model_vars(fit_lifetime, JOINT_MEDIATOR, JOINT_OUTCOME,
                 what = "The joint developmental model")

message("Joint model uses mediator: ", JOINT_MEDIATOR,
        " | outcome: ", JOINT_OUTCOME)


# ==============================================================================
# 2. Collinearity check
# ==============================================================================
# Six overlapping predictors from the same instrument. Before interpreting any
# single path, check how much they share.

period_vars <- c("TAQ_Pos_a", "TAQ_Pos_b", "TAQ_Pos_c",
                 "TAQ_Neg_a", "TAQ_Neg_b", "TAQ_Neg_c")

cor_periods <- cor(swdata[, period_vars], use = "pairwise.complete.obs")
print(round(cor_periods, 2))

vif_vals <- sapply(seq_along(period_vars), function(i) {
  r2 <- summary(lm(
    as.formula(paste(period_vars[i], "~",
                     paste(period_vars[-i], collapse = " + "))),
    data = swdata
  ))$r.squared
  1 / (1 - r2)
})
names(vif_vals) <- period_vars
print(round(vif_vals, 2))

if (any(vif_vals > 5)) {
  warning("VIF above 5 for: ",
          paste(names(vif_vals)[vif_vals > 5], collapse = ", "),
          ". Individual path estimates will be unstable.")
}


# ==============================================================================
# 3. Joint model
# ==============================================================================

# Same estimation settings as fit_mediation() in 00_functions.R.
# a1-a3 = adaptive childhood / youth / adulthood -> resilience
# b1-b3 = adverse  childhood / youth / adulthood -> resilience
# c     = resilience -> depressive symptoms (single path, shared)
# d1-d3, e1-e3 = direct paths to depressive symptoms

model_joint <- sprintf(
  "
  %1$s ~ a1*TAQ_Pos_a + a2*TAQ_Pos_b + a3*TAQ_Pos_c +
         b1*TAQ_Neg_a + b2*TAQ_Neg_b + b3*TAQ_Neg_c

  %2$s ~ c*%1$s
  %2$s ~ d1*TAQ_Pos_a + d2*TAQ_Pos_b + d3*TAQ_Pos_c +
         e1*TAQ_Neg_a + e2*TAQ_Neg_b + e3*TAQ_Neg_c

  # indirect effects per period
  ac.child := a1*c
  ac.youth := a2*c
  ac.adult := a3*c
  bc.child := b1*c
  bc.youth := b2*c
  bc.adult := b3*c

  # contrasts between periods
  ac.adult.vs.child := (a3*c) - (a1*c)
  ac.adult.vs.youth := (a3*c) - (a2*c)
  ac.youth.vs.child := (a2*c) - (a1*c)
  bc.adult.vs.child := (b3*c) - (b1*c)
  bc.adult.vs.youth := (b3*c) - (b2*c)
  bc.youth.vs.child := (b2*c) - (b1*c)
  ",
  JOINT_MEDIATOR, JOINT_OUTCOME
)

cat(model_joint)

set.seed(SEED)
fit_joint <- lavaan::sem(
  model_joint,
  data      = swdata,
  fixed.x   = FALSE,
  std.ov    = TRUE,
  missing   = "fiml",
  se        = "bootstrap",
  bootstrap = N_BOOTSTRAP
)

print(summary(fit_joint))

message("df: ", lavaan::fitMeasures(fit_joint, "df"))
print(lavaan::inspect(fit_joint, "r2"))


# ==============================================================================
# 4. Table S9
# ==============================================================================

est <- lavaan::parameterEstimates(fit_joint)

labels_joint <- c(
  ac.child = "Indirect, adaptive childhood",
  ac.youth = "Indirect, adaptive youth",
  ac.adult = "Indirect, adaptive adulthood",
  bc.child = "Indirect, adverse childhood",
  bc.youth = "Indirect, adverse youth",
  bc.adult = "Indirect, adverse adulthood",
  ac.adult.vs.child = "Contrast, adaptive: adulthood - childhood",
  ac.adult.vs.youth = "Contrast, adaptive: adulthood - youth",
  ac.youth.vs.child = "Contrast, adaptive: youth - childhood",
  bc.adult.vs.child = "Contrast, adverse: adulthood - childhood",
  bc.adult.vs.youth = "Contrast, adverse: adulthood - youth",
  bc.youth.vs.child = "Contrast, adverse: youth - childhood"
)

rows <- est[est$op == ":=", ]
rows$Term <- labels_joint[rows$lhs]

# BH correction, applied separately to the six indirect effects and to the six
# contrasts: they answer different questions.
is_contrast <- grepl("\\.vs\\.", rows$lhs)
rows$p_BH <- NA_real_
rows$p_BH[!is_contrast] <- p.adjust(rows$pvalue[!is_contrast], method = "BH")
rows$p_BH[ is_contrast] <- p.adjust(rows$pvalue[ is_contrast], method = "BH")

tableS9 <- data.frame(
  Term       = rows$Term,
  Estimate   = round(rows$est, 3),
  SE         = round(rows$se, 3),
  z          = round(rows$z, 3),
  p          = round(rows$pvalue, 3),
  `CI lower` = round(rows$ci.lower, 3),
  `CI upper` = round(rows$ci.upper, 3),
  `p (BH)`   = round(rows$p_BH, 3),
  check.names = FALSE
)

print(tableS9)

save_table(tableS9, "TableS9_joint_developmental.docx")
write.csv(tableS9, file.path(OUTPUT_DIR, "TableS9_joint_developmental.csv"),
          row.names = FALSE)


# ==============================================================================
# 5. Log
# ==============================================================================

joint_log <- c(
  "Joint developmental model with contrasts between periods",
  "All six period scores modelled simultaneously; single mediator path.",
  "Estimates are unique contributions with the other periods held constant",
  "and are therefore not comparable in size to the separate models in Table 4.",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Seed:", SEED, "| Bootstrap replications:", N_BOOTSTRAP),
  "",
  "--- Correlations among period scores ---",
  capture.output(print(round(cor_periods, 2))),
  "",
  "--- Variance inflation factors ---",
  capture.output(print(round(vif_vals, 2))),
  "",
  "--- Degrees of freedom and explained variance ---",
  paste("df:", lavaan::fitMeasures(fit_joint, "df")),
  capture.output(print(lavaan::inspect(fit_joint, "r2"))),
  "",
  "--- Table S9: indirect effects and contrasts ---",
  capture.output(print(tableS9)),
  "",
  "--- Full model summary ---",
  capture.output(summary(fit_joint)),
  "",
  "--- Session info ---",
  capture.output(sessionInfo())
)

writeLines(joint_log, file.path(OUTPUT_DIR, "joint_developmental_log.txt"))
message("Written: ", file.path(OUTPUT_DIR, "joint_developmental_log.txt"))
