# ==============================================================================
# run_all.R
#
# Runs the full analysis in the required order, in a single session.
#
# Scripts 03-05 build on data frames and fitted models left in the environment
# by 01 and 02 and cannot be run on their own; each one checks for what it
# needs and stops with a message naming the script that must run first.
#
# Run from the project root:
#   Rscript run_all.R      or      source("run_all.R")
#
# Runtime is dominated by the bootstrap (N_BOOTSTRAP = 10000 in 00_functions.R)
# and is roughly half an hour on a current desktop machine.
# ==============================================================================

scripts <- c(
  "01_main_discovery_analysis.R",        # discovery sample; Tables 1-4, S1, S2, S4
  "02_external_validation_analysis.R",   # validation sample; Table 1, S5
  "03_joint_developmental_model.R",      # needs 01; Table S9
  "04_covariate_sensitivity_discovery.R",# needs 01; Table S7
  "05_covariate_sensitivity_validation.R"# needs 02
)

for (script in scripts) {
  message("\n===== ", script, " =====")
  source(script, echo = FALSE)
}

message("\nAll scripts completed.")
