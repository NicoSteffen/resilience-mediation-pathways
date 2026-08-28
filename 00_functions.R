# ==============================================================================
# 00_functions.R
#
# Shared setup, constants and helper functions.
#
# Project : Associations between adaptive life events and depressive symptoms
#           are differentially mediated by resilience across the lifespan
#
# Sourced by every analysis script. See README.md for the run order, or use
# run_all.R.
# ==============================================================================


# ---- Packages ----------------------------------------------------------------
# Attached (used throughout):
library(tidyverse)   # data wrangling
library(lavaan)      # path / mediation models

# Called via :: (single or few uses, avoids masking):
#   foreign::read.spss     psych::alpha        broom::tidy
#   Hmisc::rcorr           report::report_sample
#   flextable::flextable   flextable::save_as_docx

options(scipen = 999)


# ---- Global analysis constants -----------------------------------------------

# Bootstrap iterations for the mediation models, as reported in the manuscript.
N_BOOTSTRAP <- 10000

# Seed for the bootstrap, applied to every model.
SEED <- 123

# Maximum number of missing CD-RISC items still tolerated (case is retained and
# person-mean imputed). The manuscript states ">10% missing" as the exclusion
# rule, which for the 25-item CD-RISC corresponds to 2.
MAX_MISSING_CDRISC <- 2


# ---- Script dependencies -----------------------------------------------------

#' Abort unless objects left in the environment by an earlier script exist
#'
#' @param ... names of the required objects
#' @param after script that must be run first, named in the error message
require_objects <- function(..., after) {
  needed  <- c(...)
  missing <- needed[!vapply(needed, exists, logical(1), envir = globalenv())]
  if (length(missing) > 0) {
    stop("Missing object(s): ", paste(missing, collapse = ", "),
         ".\n  Run ", after, " first, or run the whole pipeline with run_all.R.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Abort unless a mediator / outcome pair matches an already fitted model
#'
#' Guards against a follow-up analysis being estimated with a different
#' mediator than the main model it is compared against (e.g. the CD-RISC-10
#' of the validation sample against the CD-RISC-25 of the discovery sample).
check_model_vars <- function(fit, mediator, outcome, what = "This model") {
  vars <- lavaan::lavNames(fit, "ov")
  if (!all(c(mediator, outcome) %in% vars)) {
    stop(what, " uses mediator '", mediator, "' and outcome '", outcome,
         "', which do not both occur in the reference model (",
         paste(vars, collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}


# ---- Generic helpers ---------------------------------------------------------

#' Round every numeric column of a data frame
round_numeric <- function(df, digits = 3) {
  dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), \(x) round(x, digits)))
}

#' Format p-values with significance stars
#'
#' @param p numeric vector of p-values
#' @return character vector, e.g. "0.023*", "<.001***"
p_label <- function(p, digits = 3) {
  p <- suppressWarnings(as.numeric(p))
  stars <- dplyr::case_when(
    is.na(p)  ~ NA_character_,
    p < .001  ~ "***",
    p < .01   ~ "**",
    p < .05   ~ "*",
    TRUE      ~ ""
  )
  out <- paste0(formatC(p, format = "f", digits = digits), stars)
  out[!is.na(p) & p < .001] <- "<.001***"
  out[is.na(p)] <- NA_character_
  out
}

#' Write a data frame to a Word table in OUTPUT_DIR
#'
#' OUTPUT_DIR is set by each analysis script, so that every script writes to
#' its own folder.
save_table <- function(df, filename, dir = OUTPUT_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, filename)
  flextable::save_as_docx(flextable::flextable(as.data.frame(df)), path = path)
  message("Written: ", path)
  invisible(path)
}

#' Mean +/- SD as a formatted string
mean_sd <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"),
          mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

#' n (%) for each level of a factor / character vector
n_percent <- function(x, digits = 1) {
  tab <- table(x, useNA = "no")
  out <- sprintf(paste0("%d (%.", digits, "f%%)"), tab, 100 * tab / sum(tab))
  stats::setNames(out, names(tab))
}


# ---- TAQ scoring -------------------------------------------------------------

# Scoring follows the TAQ manual (Spinazzola et al., 2019). The German version
# administered here excludes the two "family secrets" items and the check-up
# item 43, leaving 40 items renumbered 1-40, over three developmental periods
# (a = childhood 7-12, b = youth 13-18, c = adulthood >18).
#
# Response-label -> score maps live in each analysis script, because the two
# samples were collected with slightly different SPSS value labels. The scoring
# itself is identical:
#   often / very much = 3, occasionally = 2, rarely = 0, never = 0,
#   yes = 3, no = 0, don't know = 0
# and reversed for TAQ item 2, the only positively worded item.

#' Recode item columns from labelled responses to numeric scores
#'
#' Unmatched labels become NA and are reported, so that silent scoring
#' failures cannot pass unnoticed.
recode_items <- function(data, cols, map) {
  observed <- unique(unlist(lapply(data[cols], \(x) as.character(x))))
  observed <- observed[!is.na(observed)]
  unmatched <- setdiff(observed, names(map))
  if (length(unmatched) > 0) {
    warning("Unmatched response labels set to NA: ",
            paste(sQuote(unmatched), collapse = ", "), call. = FALSE)
  }
  data[cols] <- lapply(data[cols], \(x) unname(map[as.character(x)]))
  data
}

#' TAQ subscale item assignment, by developmental period
#'
#' a = childhood (7-12), b = youth (13-18), c = adulthood (>18)
TAQ_SUBSCALES <- local({
  items <- list(
    TAQ_Competence        = c(3, 4),
    TAQ_Safety            = c(1, 5, 8),
    TAQ_Neglect           = c(2, 6, 7, 21, 27),
    TAQ_Separation        = c(10, 11, 12, 14),
    TAQ_Emotional_Abuse   = c(9, 16, 17, 18, 19),
    TAQ_Physical_Abuse    = c(28, 29, 30),
    TAQ_Sexual_Abuse      = c(35, 36, 37, 38),
    TAQ_Witnessing        = c(20, 22, 23, 24, 31, 34),
    TAQ_Impersonal_Trauma = c(13, 15, 32, 33, 39, 40),
    TAQ_Alcohol_Drugs     = c(25, 26)
  )
  lapply(items, \(nums) {
    stats::setNames(
      lapply(c("a", "b", "c"), \(p) paste0("TAQ_", nums, p)),
      c("a", "b", "c")
    )
  })
})

# Adaptive (TAQ+) vs. adverse (TAQ-) subscales
TAQ_POSITIVE <- c("TAQ_Safety", "TAQ_Competence")
TAQ_NEGATIVE <- c("TAQ_Neglect", "TAQ_Separation", "TAQ_Emotional_Abuse",
                  "TAQ_Physical_Abuse", "TAQ_Sexual_Abuse", "TAQ_Witnessing",
                  "TAQ_Impersonal_Trauma", "TAQ_Alcohol_Drugs")

#' Subscale score for one developmental period
#'
#' Mean of the ENDORSED items only: a response of 0 is treated as
#' "not endorsed" rather than as a score of zero. A respondent who endorsed
#' nothing scores 0; a respondent with no data at all scores NaN.
calc_subscale <- function(data, items) {
  apply(data[, items, drop = FALSE], 1, function(row) {
    if (all(is.na(row))) {
      NaN
    } else if (all(is.na(row) | row == 0)) {
      0
    } else {
      mean(ifelse(row == 0, NA, row), na.rm = TRUE)
    }
  })
}

#' Add all TAQ subscale scores, period totals and TAQ+ / TAQ- composites
#'
#' Creates, for each subscale: <name>_a, <name>_b, <name>_c and <name>_total,
#' plus TAQ_Pos / TAQ_Neg (lifetime, range 0-9) and
#' TAQ_Pos_a/b/c, TAQ_Neg_a/b/c (period-specific, range 0-3).
#'
#' Period totals are deliberately NOT computed with na.rm = TRUE: if a
#' respondent did not complete the items for a developmental period, that
#' period must stay missing rather than count as an absence of events, which
#' would produce an artificially low lifetime score. Cases are retained;
#' missingness at the composite level is handled by FIML in the path models.
add_taq_scores <- function(data) {
  for (subscale in names(TAQ_SUBSCALES)) {
    for (period in c("a", "b", "c")) {
      data[[paste0(subscale, "_", period)]] <-
        calc_subscale(data, TAQ_SUBSCALES[[subscale]][[period]])
    }
    data[[paste0(subscale, "_total")]] <-
      rowSums(data[, paste0(subscale, "_", c("a", "b", "c"))])
  }

  # Lifetime composites (sum across periods, then averaged over subscales)
  data$TAQ_Pos <- rowMeans(data[, paste0(TAQ_POSITIVE, "_total")])
  data$TAQ_Neg <- rowMeans(data[, paste0(TAQ_NEGATIVE, "_total")])

  # Period-specific composites
  for (period in c("a", "b", "c")) {
    data[[paste0("TAQ_Pos_", period)]] <-
      rowMeans(data[, paste0(TAQ_POSITIVE, "_", period)], na.rm = TRUE)
    data[[paste0("TAQ_Neg_", period)]] <-
      rowMeans(data[, paste0(TAQ_NEGATIVE, "_", period)], na.rm = TRUE)
  }

  n_incomplete <- sum(is.na(data$TAQ_Pos))
  if (n_incomplete > 0) {
    message("Lifetime TAQ composites set to NA for ", n_incomplete,
            " case(s) with an incomplete developmental period.")
  }
  data
}


# ---- CD-RISC scoring ---------------------------------------------------------

#' Person-mean imputation of CD-RISC items
#'
#' Applied only to cases retained by the MAX_MISSING_CDRISC rule.
impute_person_mean <- function(data, cols) {
  mat <- as.matrix(data[, cols, drop = FALSE])
  row_means <- rowMeans(mat, na.rm = TRUE)
  na_idx <- which(is.na(mat), arr.ind = TRUE)
  mat[na_idx] <- row_means[na_idx[, "row"]]
  data[, cols] <- as.data.frame(mat)
  data
}


# ---- Mediation models --------------------------------------------------------

#' Fit the resilience mediation model
#'
#' Adaptive (x_pos) and adverse (x_neg) life events -> mediator -> outcome,
#' with direct paths retained. Used for the lifetime, period-specific and
#' validation models alike.
#'
#' Mediator and outcome are arguments, never taken from the global
#' environment, so that a model cannot silently inherit the wrong scale.
fit_mediation <- function(data, x_pos, x_neg, mediator, outcome,
                          seed = SEED, n_boot = N_BOOTSTRAP) {
  model <- sprintf(
    "
    %1$s ~ a*%2$s + b*%3$s
    %4$s ~ c*%1$s
    %4$s ~ d*%2$s + e*%3$s
    ac        := a*c
    bc        := b*c
    total.pos := d + (a*c)
    total.neg := e + (b*c)
    ",
    mediator, x_pos, x_neg, outcome
  )

  set.seed(seed)
  lavaan::sem(
    model,
    data      = data,
    fixed.x   = FALSE,
    std.ov    = TRUE,
    missing   = "fiml",
    se        = "bootstrap",
    bootstrap = n_boot
  )
}

#' Tidy a fitted mediation model into a publication table
#'
#' Keeps regression paths, the covariance between the two exogenous
#' predictors, and the defined (indirect / total) effects. Drops intercepts
#' and all variance terms.
#'
#' @param labels named character vector, variable name -> display label.
#'   Order matters: replacement is sequential, so list longer names first
#'   (e.g. "TAQ_Pos_a" before "TAQ_Pos").
tidy_mediation <- function(fit, labels, digits = 3) {
  est <- as.data.frame(lavaan::parameterEstimates(fit))

  defined <- c(ac        = "TAQ+ indirect",
               bc        = "TAQ- indirect",
               total.pos = "TAQ+ total",
               total.neg = "TAQ- total")

  term <- ifelse(
    est$op == ":=",
    defined[est$lhs],
    paste(est$lhs, ifelse(est$op == "~~", "\u2194", "\u2190"), est$rhs)
  )
  for (nm in names(labels)) {
    term <- gsub(nm, labels[[nm]], term, fixed = TRUE)
  }

  keep <- (est$op %in% c("~", ":=")) |
    (est$op == "~~" & est$lhs != est$rhs)   # covariance, not variance

  out <- est[keep, c("est", "se", "z", "pvalue", "ci.lower", "ci.upper")]
  out <- tibble::add_column(out, Term = term[keep], .before = 1)
  out <- round_numeric(out, digits)

  names(out) <- c("Term", "Point estimate", "Standard error", "z",
                  "p-value", "CI lower", "CI upper")
  out$`p-value` <- p_label(out$`p-value`, digits)
  rownames(out) <- NULL
  out
}


# ---- Reliability -------------------------------------------------------------

#' Cronbach's alpha with mean inter-item correlation
#'
#' Items with zero variance, and items whose pairwise correlations remain
#' undefined because of near-zero endorsement, are dropped and listed.
#' Returns NA where fewer than two usable items remain.
alpha_summary <- function(items, label, period = NA_character_) {
  items <- items[, colSums(!is.na(items)) > 0, drop = FALSE]

  # 1. zero-variance items
  v <- vapply(items, var, numeric(1), na.rm = TRUE)
  dropped <- names(items)[is.na(v) | v == 0]
  items <- items[, !names(items) %in% dropped, drop = FALSE]

  # 2. items still producing undefined correlations - drop one at a time,
  #    worst first, since undefined cells are pairwise
  while (ncol(items) >= 2) {
    R <- suppressWarnings(cor(items, use = "pairwise.complete.obs"))
    na_counts <- colSums(is.na(R))
    if (all(na_counts == 0)) break
    worst <- names(which.max(na_counts))
    dropped <- c(dropped, worst)
    items <- items[, names(items) != worst, drop = FALSE]
  }

  k <- ncol(items)
  if (k < 2) {
    a <- NA_real_
    r_bar <- NA_real_
  } else {
    a <- suppressWarnings(
      psych::alpha(items, use = "pairwise", warnings = FALSE)$total$raw_alpha
    )
    r_bar <- a / (k - a * (k - 1))
  }

  tibble::tibble(
    Scale               = label,
    Period              = period,
    k                   = k,
    Alpha               = round(a, 2),
    `Mean inter-item r` = round(r_bar, 3),
    `Items dropped`     = if (length(dropped)) paste(dropped, collapse = ", ") else "-"
  )
}

#' Item columns for a TAQ composite
#'
#' @param which "positive", "negative" or a subscale name from TAQ_SUBSCALES
#' @param period "a", "b", "c", or "all" for the pooled item set
taq_items <- function(data, which, period = "all") {
  subscales <- switch(which,
                      positive = TAQ_POSITIVE,
                      negative = TAQ_NEGATIVE,
                      which
  )
  periods <- if (identical(period, "all")) c("a", "b", "c") else period
  cols <- unlist(lapply(subscales, \(s) unlist(TAQ_SUBSCALES[[s]][periods])))
  data[, intersect(cols, names(data)), drop = FALSE]
}

#' Full reliability table for one sample
#'
#' Reports Cronbach's alpha and the mean inter-item correlation for the TAQ
#' composites (pooled and by developmental period), for each TAQ subscale by
#' period, and for the CD-RISC and BDI-II.
reliability_table <- function(data, cdrisc_cols, bdi_cols) {
  periods <- c(a = "Childhood", b = "Youth", c = "Adulthood")

  rows <- list(
    alpha_summary(taq_items(data, "positive"), "TAQ+ (adaptive)", "Lifetime"),
    alpha_summary(taq_items(data, "negative"), "TAQ- (adverse)",  "Lifetime")
  )

  for (p in names(periods)) {
    rows <- c(rows, list(
      alpha_summary(taq_items(data, "positive", p), "TAQ+ (adaptive)", periods[[p]]),
      alpha_summary(taq_items(data, "negative", p), "TAQ- (adverse)",  periods[[p]])
    ))
  }

  for (s in names(TAQ_SUBSCALES)) {
    label <- sub("^TAQ_", "", s)
    label <- gsub("_", " ", label)
    rows <- c(rows, list(alpha_summary(taq_items(data, s), label, "Lifetime")))
    for (p in names(periods)) {
      rows <- c(rows, list(alpha_summary(taq_items(data, s, p), label, periods[[p]])))
    }
  }

  rows <- c(rows, list(
    alpha_summary(data[, cdrisc_cols, drop = FALSE], "CD-RISC", "-"),
    alpha_summary(data[, bdi_cols,    drop = FALSE], "BDI-II",  "-")
  ))

  dplyr::bind_rows(rows)
}
