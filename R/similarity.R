#' Multi-group tests of profile similarity
#'
#' Fits the sequence of configural, structural, dispersion and distributional
#' similarity models of Morin et al. (2016) for a given number of profiles and
#' applies the decision rule used in the MG-LPA Analyzer: similarity at a
#' given level is supported when at least two of the three information
#' criteria (CAIC, BIC and ABIC) decrease relative to the previous, less
#' constrained model; it is rejected when none decreases; otherwise the
#' evidence is mixed.
#'
#' @inheritParams enumerate_profiles
#' @param group Name of the known-group variable.
#' @param k Number of profiles retained in every group at the enumeration
#'   stage.
#' @param models Which models to estimate, in order.
#' @return An object of class `mglpa_similarity` with elements `fits` (named
#'   list of `mglpa_fit` objects), `table` (fit indices, differences and
#'   decisions), `retained` (name of the most constrained model supported by
#'   the sequence) and `data`.
#' @examples
#' \donttest{
#' data(mglpa_wm)
#' sim <- similarity_tests(mglpa_wm,
#'   indicators = c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed"),
#'   group = "group", k = 3, starts = c(20, 4), seed = 1)
#' sim
#' }
#' @export
similarity_tests <- function(data, indicators, group, k,
                             variances = c("class_invariant", "class_varying"),
                             models = c("configural", "structural", "dispersion", "distributional"),
                             starts = c(100, 10), stiterations = 10, maxiter = 500, tol = 1e-6,
                             se = c("robust", "ml", "none"),
                             order_profiles = c("level", "size", "none"),
                             group_levels = NULL, seed = NULL, verbose = TRUE, cores = 1L) {
  variances <- match.arg(variances)
  se <- match.arg(se)
  order_profiles <- match.arg(order_profiles)
  models <- match.arg(models, several.ok = TRUE)
  if (!is.null(seed)) set.seed(seed)
  specs <- list(configural = c(FALSE, FALSE, FALSE), structural = c(TRUE, FALSE, FALSE),
                dispersion = c(TRUE, TRUE, FALSE), distributional = c(TRUE, TRUE, TRUE))
  fits <- list()
  prev <- NULL
  for (nm in models) {
    cc <- specs[[nm]]
    model <- mglpa_model(k, indicators, group = group, variances = variances,
                         equal_means = cc[1], equal_variances = cc[2],
                         equal_proportions = cc[3], group_levels = group_levels)
    f <- fit_lpa(data, model, starts = starts, stiterations = stiterations,
                 maxiter = maxiter, tol = tol, start_values = prev, se = se, verbose = FALSE,
                 cores = cores)
    if (order_profiles != "none") f <- reorder_profiles(f, order_profiles)
    fits[[nm]] <- f
    prev <- f
    if (verbose)
      message(sprintf("%-14s LL = %.3f, npar = %d, BIC = %.3f, CAIC = %.3f, ABIC = %.3f, best LL replicated %d/%d",
                      nm, f$loglik, f$npar, f$fit[["BIC"]], f$fit[["CAIC"]], f$fit[["ABIC"]],
                      f$starts$best_ll_replications, length(f$starts$final_ll)))
  }
  table <- similarity_table_from_fits(fits)
  retained <- retained_model(table)
  structure(list(fits = fits, table = table, retained = retained, indicators = indicators,
                 group = group, k = k, variances = variances, data = data, call = match.call()),
            class = "mglpa_similarity")
}

#' Compare two nested models with the two-out-of-three rule
#'
#' @param new,old `mglpa_fit` objects; `new` is the more constrained model.
#' @return A list with the differences in CAIC, BIC and ABIC (`new - old`),
#'   the number of criteria that decreased and the decision (`"supported"`,
#'   `"mixed"` or `"rejected"`).
#' @export
compare_fit <- function(new, old) {
  ics <- c("CAIC", "BIC", "ABIC")
  diff <- new$fit[ics] - old$fit[ics]
  n <- sum(diff < 0)
  decision <- if (n >= 2) "supported" else if (n == 0) "rejected" else "mixed"
  list(diff = diff, n_decrease = n, decision = decision)
}

similarity_table_from_fits <- function(fits) {
  rows <- lapply(names(fits), function(nm) fit_row(fits[[nm]], model = nm))
  tab <- do.call(rbind, rows)
  tab$group <- NULL
  tab$dCAIC <- tab$dBIC <- tab$dABIC <- NA_real_
  tab$n_decrease <- NA_integer_
  tab$decision <- NA_character_
  if (length(fits) > 1L) for (i in 2:length(fits)) {
    cmp <- compare_fit(fits[[i]], fits[[i - 1L]])
    tab$dCAIC[i] <- cmp$diff[["CAIC"]]
    tab$dBIC[i] <- cmp$diff[["BIC"]]
    tab$dABIC[i] <- cmp$diff[["ABIC"]]
    tab$n_decrease[i] <- cmp$n_decrease
    tab$decision[i] <- cmp$decision
  }
  rownames(tab) <- NULL
  tab
}

retained_model <- function(tab) {
  retained <- tab$model[1L]
  if (nrow(tab) > 1L) for (i in 2:nrow(tab)) {
    if (identical(tab$decision[i], "supported")) retained <- tab$model[i] else break
  }
  retained
}

#' @export
print.mglpa_similarity <- function(x, digits = 3, ...) {
  cat("<mglpa_similarity> ", x$k, "-profile solution, ", length(x$fits), " model(s), ",
      gsub("_", "-", x$variances), " variances\n", sep = "")
  tab <- x$table
  out <- data.frame(model = tab$model, LL = round(tab$LL, digits), npar = tab$npar,
                    AIC = round(tab$AIC, digits), BIC = round(tab$BIC, digits),
                    CAIC = round(tab$CAIC, digits), ABIC = round(tab$ABIC, digits),
                    entropy = round(tab$entropy, digits),
                    decrease = ifelse(is.na(tab$n_decrease), "", paste0(tab$n_decrease, "/3")),
                    decision = ifelse(is.na(tab$decision), "", tab$decision),
                    reps = paste0(tab$replications, "/", tab$n_final),
                    conv = ifelse(tab$converged, "yes", "NO"))
  print(out, row.names = FALSE)
  cat("\nDecision rule: at least 2 of CAIC, BIC and ABIC lower than the previous model.\n")
  cat("Most constrained model supported: ", x$retained, "\n", sep = "")
  cat("\nProfile proportions by group (", x$retained, " model):\n", sep = "")
  print(round(x$fits[[x$retained]]$class_counts$by_group_model, digits))
  print_estimation_notes(x$fits)
  invisible(x)
}

#' Fit-index table of the similarity tests
#'
#' @param x An `mglpa_similarity` object.
#' @return A data frame.
#' @export
similarity_table <- function(x) {
  stopifnot(inherits(x, "mglpa_similarity"))
  x$table
}
