# ==============================================================================
# 02_external_validation_analysis.R
#
# Independent validation sample: N = 67 patients with persistent depressive
# disorder (PDD, coded "CD") or borderline personality disorder (BPD),
# collected at the same department for a separate study (Barton et al., 2023).
# Healthy controls from that study are excluded here.
#
# Applies the lifetime mediation model of 01_main_discovery_analysis.R with
# identical specification. Resilience is measured with the CD-RISC-10 rather
# than the CD-RISC-25.
#
# Produces, in outputs_external_validation/:
#   Reliability_validation.docx        internal consistency
#   Table1_sample_validation.docx      sample characteristics
#   TableS5_mediation_validation.docx  validation mediation model
#
# Independent of 01; requires nothing but 00_functions.R.
# Run from the project root. See README.md for the run order.
# ==============================================================================

source("00_functions.R")

OUTPUT_DIR <- "outputs_external_validation"

# SPSS export of the validation sample. Not published (patient data);
# see README.md for access conditions.
DATA_FILE <- file.path("raw_data", "file.sav")

# Grouping variable identifying the two patient groups. The sample size is
# asserted below (PDD n = 33, BPD n = 34).
GROUP_VAR <- "group"


# ==============================================================================
# 1. Import
# ==============================================================================

val <- foreign::read.spss(DATA_FILE, to.data.frame = TRUE)


# ==============================================================================
# 2. Preprocessing
# ==============================================================================

## 2.1 Drop instruments not used in this paper ---------------------------------

drop_prefixes <- c("KOMME", "NEO_FFI", "IPDE_", "MADRS_", "BSL_", "IIPC_",
                   "RSQ", "BIS_", "UCLA_", "SNI_", "BUL_", "STAI_", "VAR",
                   "X.", "ES", "FNS", "Behav", "NTS", "pre_", "post_", "p2",
                   "soe", "d_", "SICD")

drop_cols <- c(
  paste0("CDRISK_", 1:25, "_alterFB"),   # age-of-event version, not used
  "TAQ_43",                              # check-up item, not part of any subscale
  unlist(lapply(drop_prefixes, \(p) grep(paste0("^", p), names(val), value = TRUE)))
)

message("Dropping ", length(intersect(drop_cols, names(val))), " unused columns.")
val <- val[, !names(val) %in% drop_cols]


## 2.2 CD-RISC completeness ----------------------------------------------------

cdrisc_cols <- grep("^CDRISK_", names(val), value = TRUE)

# The validation sample used the 10-item CD-RISC.
stopifnot(length(cdrisc_cols) == 10)

missing_cdrisc <- rowSums(is.na(val[, cdrisc_cols]))
message("CD-RISC missing items per case:")
print(table(missing_cdrisc))

val <- val[missing_cdrisc <= MAX_MISSING_CDRISC, ]


## 2.3 Analysis sample: PDD ("CD") and BPD only --------------------------------

message("Group sizes before selection:")
print(table(val[[GROUP_VAR]]))

val <- val[val[[GROUP_VAR]] %in% c("BPD", "CD"), ]
val[[GROUP_VAR]] <- droplevels(val[[GROUP_VAR]])

message("Analysis sample: n = ", nrow(val))
stopifnot(nrow(val) == 67)


## 2.4 TAQ item selection and renumbering --------------------------------------
# As in the discovery sample: items 9 and 24 ("family secrets") are excluded
# and the remaining 40 items x 3 developmental periods renumbered 1a-40c.

val <- val[, !grepl("TAQ_9|TAQ_24", names(val))]

taq_idx <- grep("TAQ", names(val))
stopifnot(length(taq_idx) == 120)
names(val)[taq_idx] <- paste0("TAQ_", rep(1:40, each = 3),
                              rep(letters[1:3], times = 40))


# ==============================================================================
# 3. Medication coding
# ==============================================================================
# Each substance is a separate dose column; a non-missing dose means the drug
# was prescribed at admission.

antidepressants <- c(
  "Agomelantin",        # agomelatine
  "Amitriptylin",       # tricyclic
  "Amioxid",            # amoxapine
  "Bupropion_Dosis", "Citalopram_Dosis", "Doxepin_Dosis", "Duloxetin",
  "Escitalopram_Dosis",
  "Fluoextin_Dosis",    # fluoxetine (spelling as in the export)
  "Mitrazapin_Dosis",   # mirtazapine (spelling as in the export)
  "Milnacipran", "Paroxetin_Dosis", "Sertralin",
  "Tianeurax",          # tianeptine
  "Trazodon", "Trimipramin", "Tranylcypromin"
)

antipsychotics <- c(
  "Aripiprazol", "Melperon_Dosis", "Olanzapin", "Quetiapin_Dosis",
  "Risperidon", "Promethazin", "Prothipendyl"
)

mood_stabilizers <- c(
  "Ergenyl",        # valproate
  "Lamictal",       # lamotrigine
  "Lithium", "Levetiracetam", "Oxcarbazepin", "Gabapentin", "Pregabalin"
)

other_psychotropics <- c(
  "adderall", "Methylphenidat", "Modafinil", "Clonidin",
  "Circadin",       # melatonin
  "Naltrexon"
)

psychotropics <- c(antidepressants, antipsychotics,
                   mood_stabilizers, other_psychotropics)

# Guard against silent misses if a substance column is absent from the export.
missing_cols <- setdiff(psychotropics, names(val))
if (length(missing_cols) > 0) {
  warning("Medication columns not found and ignored: ",
          paste(missing_cols, collapse = ", "), call. = FALSE)
}

any_prescribed <- function(data, cols) {
  cols <- intersect(cols, names(data))
  ifelse(rowSums(!is.na(data[, cols, drop = FALSE])) > 0, "yes", "no")
}

val$medication      <- any_prescribed(val, psychotropics)
val$antidepressants <- any_prescribed(val, antidepressants)
val$antipsychotics  <- any_prescribed(val, antipsychotics)


# ==============================================================================
# 4. TAQ scoring
# ==============================================================================
# Same scoring as the discovery sample; the SPSS value labels differ slightly.
#
# Nine yes/no TAQ columns carry the undeclared numeric codes "5" and "3"
# instead of a value label. Per the value-label table of this export they
# follow the Likert scheme of the remaining TAQ items (1 = never, 2 = rarely,
# 3 = occasionally, 4 = often, 5 = don't know), so they are scored as
# "occasionally" and "don't know" rather than taken at face value.

TAQ_MAP <- c(
  "oft/stark"                    = 3,
  "Gelegentlich/im Grunde schon" = 2,
  "selten/ein wenig"             = 0,
  "nie/gar nicht"                = 0,
  "ja"                           = 3,
  "nein"                         = 0,
  "weiß nicht"                   = 0,
  "5"                            = 0,   # don't know
  "3"                            = 2    # occasionally
)

TAQ_MAP_REVERSED <- c(
  "oft/stark"                    = 0,
  "Gelegentlich/im Grunde schon" = 0,
  "selten/ein wenig"             = 2,
  "nie/gar nicht"                = 3,
  "weiß nicht"                   = 0,
  "5"                            = 0,   # don't know
  "3"                            = 0    # occasionally, reverse-scored
)

# TAQ item 2 is the only reverse-scored item; the anchored pattern keeps
# items 20-29 out. See README.md, "Reverse-scoring of TAQ item 2".
reverse_cols <- grep("^TAQ_2[abc]$", names(val), value = TRUE)
message("Reverse-scored TAQ columns (n = ", length(reverse_cols), "): ",
        paste(reverse_cols, collapse = ", "))

standard_cols <- setdiff(grep("^TAQ_[0-9]", names(val), value = TRUE),
                         reverse_cols)

# Reverse-scored columns first: recoding turns labels into numbers, so a column
# recoded once can no longer be matched by the second map.
val <- recode_items(val, reverse_cols,  TAQ_MAP_REVERSED)
val <- recode_items(val, standard_cols, TAQ_MAP)

val <- add_taq_scores(val)


# ==============================================================================
# 5. CD-RISC-10 scoring
# ==============================================================================

CDRISC_MAP <- c(
  "überhaupt nicht wahr" = 0,
  "selten wahr"          = 1,
  "manchmal wahr"        = 2,
  "oft wahr"             = 3,
  "fast immer wahr"      = 4
)

val <- recode_items(val, cdrisc_cols, CDRISC_MAP)
val <- impute_person_mean(val, cdrisc_cols)

val$CDRISC_10_sum <- rowSums(val[, cdrisc_cols])


# ==============================================================================
# 6. Completeness of the TAQ composites and internal consistency
# ==============================================================================

taq_n <- data.frame(
  Score = c("TAQ_Pos / TAQ_Neg (lifetime)",
            "childhood", "youth", "adulthood"),
  n = c(sum(!is.na(val$TAQ_Pos)),
        sum(!is.na(val$TAQ_Pos_a)),
        sum(!is.na(val$TAQ_Pos_b)),
        sum(!is.na(val$TAQ_Pos_c)))
)
print(taq_n)

items_bdi <- grep("^BDI_[0-9]", names(val), value = TRUE)
reliability_val <- reliability_table(val, cdrisc_cols, items_bdi)
print(reliability_val, n = Inf)
save_table(reliability_val, "Reliability_validation.docx")


# ==============================================================================
# 7. Table 1 -- validation sample characteristics
# ==============================================================================

table1_validation <- tibble::tibble(
  Variable = c(
    paste("Sex:",             names(n_percent(val$Geschlecht))),
    paste("Group:",           names(n_percent(val[[GROUP_VAR]]))),
    "Age (years)",
    "Age at onset (years)",
    paste("Medication:",      names(n_percent(val$medication))),
    paste("Antidepressants:", names(n_percent(val$antidepressants))),
    paste("Antipsychotics:",  names(n_percent(val$antipsychotics))),
    "TAQ+ (adaptive)",
    "TAQ- (adverse)",
    "CD-RISC-10",
    "BDI-II"
  ),
  Value = c(
    n_percent(val$Geschlecht),
    n_percent(val[[GROUP_VAR]]),
    mean_sd(val$Alter_echt),
    mean_sd(val$Alter_Krankheitsbeginn),
    n_percent(val$medication),
    n_percent(val$antidepressants),
    n_percent(val$antipsychotics),
    mean_sd(val$TAQ_Pos),
    mean_sd(val$TAQ_Neg),
    mean_sd(val$CDRISC_10_sum),
    mean_sd(val$BDI_Score)
  )
)

print(table1_validation, n = Inf)
save_table(table1_validation, "Table1_sample_validation.docx")


# ==============================================================================
# 8. Table S5 -- validation mediation model
# ==============================================================================

# Mediator and outcome of this sample. Note that the mediator is the
# CD-RISC-10, not the CD-RISC-25 used in the discovery sample.
VAL_MEDIATOR <- "CDRISC_10_sum"
VAL_OUTCOME  <- "BDI_Score"

stopifnot(all(c(VAL_MEDIATOR, VAL_OUTCOME) %in% names(val)))

fit_validation <- fit_mediation(
  data     = val,
  x_pos    = "TAQ_Pos",
  x_neg    = "TAQ_Neg",
  mediator = VAL_MEDIATOR,
  outcome  = VAL_OUTCOME
)
print(summary(fit_validation))

tableS5 <- tidy_mediation(
  fit_validation,
  labels = c(CDRISC_10_sum = "CDRISC", BDI_Score = "BDI",
             TAQ_Pos = "TAQ+", TAQ_Neg = "TAQ-")
)
print(tableS5)
save_table(tableS5, "TableS5_mediation_validation.docx")


# ==============================================================================
# 9. Session information
# ==============================================================================

sessionInfo()
