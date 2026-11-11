#' Explanatory similarity: profile-specific levels of distal outcomes across groups
#'
#' Tests whether the profile-specific means of distal outcomes are the same
#' across groups (step 6 of Morin et al., 2016). The profile solution,
#' including the profile-membership logits, is held fixed at the estimates
#' of a previously retained similarity model. For each residual-variance
#' specification of the outcomes (one variance per outcome, or profile- and
#' group-specific variances) a model with group-specific outcome means is
#' compared with a model whose outcome means are equal across groups.
#' Explanatory similarity is supported when at least two of CAIC, BIC and
#' ABIC are lower for the constrained model. Pairwise differences between
#' profiles (within groups) and between groups (within profiles) are obtained
#' with the delta method.
#'
#' @inheritParams predictive_similarity
#' @param outcomes Character vector of outcome names.
#' @param outcome_variances `"both"`, `"equal"` (one residual variance per
#'   outcome) or `"free"` (profile- and group-specific variances).
#' @return An object of class `mglpa_explanatory` with elements `fits`
#'   (named `free_means_equal_var`, `equal_means_equal_var`,
#'   `free_means_free_var`, `equal_means_free_var`), `table`, `decisions`,
#'   `means` (outcome means by profile and group, with confidence limits) and
#'   `contrasts` (pairwise differences), the last two by model.
#' @export
explanatory_similarity <- function(fit, outcomes, data = NULL, source = NULL,
                                   outcome_variances = c("both", "equal", "free"),
                                   starts = c(20, 5), stiterations = 10, maxiter = 500,
                                   tol = 1e-6, se = c("robust", "ml", "none"), digits = 5,
                                   seed = NULL, verbose = TRUE) {
  outcome_variances <- match.arg(outcome_variances)
  se <- match.arg(se)
  if (!is.null(seed)) set.seed(seed)
  src <- resolve_source(fit, source)
  if (is.null(data)) data <- src$data$data_frame
  if (is.null(data)) stop("`data` must be supplied.")
  src_fit <- src$fit
  m0 <- src_fit$model
  if (is.null(m0$group)) stop("Explanatory similarity requires a model with a group variable.")
  fixed <- fixed_from_fit(src_fit, c("mean", "variance", "intercept", "slope"), digits = digits)
  var_specs <- if (outcome_variances == "both") c("equal", "free") else outcome_variances
  fits <- list()
  for (ov in var_specs) for (om in c("free", "equal")) {
    model <- mglpa_model(m0$k, m0$indicators, group = m0$group, variances = m0$variances,
                         equal_means = m0$equal_means, equal_variances = m0$equal_variances,
                         equal_proportions = m0$equal_proportions, predictors = m0$predictors,
                         predictor_effects = m0$predictor_effects, outcomes = outcomes,
                         outcome_means = om, outcome_variances = ov,
                         group_levels = m0$group_levels)
    nm <- sprintf("%s_means_%s_var", om, ov)
    f <- fit_lpa(data, model, starts = starts, stiterations = stiterations, maxiter = maxiter,
                 tol = tol, fixed = fixed, se = se, verbose = FALSE)
    fits[[nm]] <- f
    if (verbose)
      message(sprintf("%-22s LL = %.3f, npar = %d, BIC = %.3f, CAIC = %.3f, ABIC = %.3f",
                      nm, f$loglik, f$npar, f$fit[["BIC"]], f$fit[["CAIC"]], f$fit[["ABIC"]]))
  }
  rows <- lapply(names(fits), function(nm) fit_row(fits[[nm]], model = nm))
  tab <- do.call(rbind, rows)
  tab$group <- NULL
  tab$decision <- NA_character_
  tab$n_decrease <- NA_integer_
  decisions <- list()
  for (ov in var_specs) {
    free_nm <- sprintf("free_means_%s_var", ov)
    eq_nm <- sprintf("equal_means_%s_var", ov)
    cmp <- compare_fit(fits[[eq_nm]], fits[[free_nm]])
    decisions[[ov]] <- cmp
    i <- which(tab$model == eq_nm)
    tab$decision[i] <- cmp$decision
    tab$n_decrease[i] <- cmp$n_decrease
  }
  means <- lapply(fits, function(f) {
    pm <- profile_means(f, standardized = FALSE)
    pm[pm$variable %in% outcomes, ]
  })
  contrasts <- lapply(fits, outcome_contrasts)
  structure(list(fits = fits, table = tab, decisions = decisions, means = means,
                 contrasts = contrasts, outcomes = outcomes, source = src$name,
                 source_fit = src_fit, data = data, call = match.call()),
            class = "mglpa_explanatory")
}

#' Pairwise differences in outcome means between profiles and between groups
#'
#' @param fit An `mglpa_fit` with distal outcomes.
#' @return A data frame of delta-method contrasts: differences between pairs
#'   of profiles within each group (or `"all"` when outcome means are equal
#'   across groups) and, when outcome means are group-specific, differences
#'   between pairs of groups within each profile.
#' @export
outcome_contrasts <- function(fit) {
  stopifnot(inherits(fit, "mglpa_fit"))
  m <- fit$model
  if (!length(m$outcomes)) stop("The model has no outcomes.")
  K <- m$k
  G <- fit$D$G
  lev <- as.character(fit$data$group_levels)
  equal <- m$outcome_means == "equal"
  lab <- function(z, k, g) if (equal) sprintf("m_%s_c%d", z, k) else sprintf("m_%s_c%d_g%d", z, k, g)
  L <- list(); meta <- list()
  groups <- if (equal) list(list(g = NA, label = "all")) else
    lapply(seq_len(G), function(g) list(g = g, label = lev[g]))
  for (z in m$outcomes) for (gr in groups) for (k1 in seq_len(K - 1L)) for (k2 in (k1 + 1L):K) {
    w <- stats::setNames(c(1, -1), c(lab(z, k1, gr$g), lab(z, k2, gr$g)))
    id <- sprintf("%s|%s|p%d-p%d", z, gr$label, k1, k2)
    L[[id]] <- w
    meta[[id]] <- data.frame(outcome = z, type = "between profiles", group_level = gr$label,
                             profile = k1, comparison = sprintf("Profile %d - Profile %d", k1, k2),
                             stringsAsFactors = FALSE)
  }
  if (!equal && G > 1L) {
    for (z in m$outcomes) for (k in seq_len(K)) for (g1 in seq_len(G - 1L)) for (g2 in (g1 + 1L):G) {
      w <- stats::setNames(c(1, -1), c(lab(z, k, g1), lab(z, k, g2)))
      id <- sprintf("%s|p%d|g%d-g%d", z, k, g1, g2)
      L[[id]] <- w
      meta[[id]] <- data.frame(outcome = z, type = "between groups", group_level = "",
                               profile = k, comparison = sprintf("Profile %d: group %s - group %s",
                                                                 k, lev[g1], lev[g2]),
                               stringsAsFactors = FALSE)
    }
  }
  ce <- contrast_estimates(fit, L)
  out <- cbind(do.call(rbind, meta), ce[, c("estimate", "se", "z", "p", "ci_lower", "ci_upper")])
  rownames(out) <- NULL
  out
}

#' @export
print.mglpa_explanatory <- function(x, digits = 3, ...) {
  cat("<mglpa_explanatory> outcomes: ", paste(x$outcomes, collapse = ", "),
      "; profile solution fixed at the ", x$source, " model\n", sep = "")
  tab <- x$table
  out <- data.frame(model = tab$model, LL = round(tab$LL, digits), npar = tab$npar,
                    AIC = round(tab$AIC, digits), BIC = round(tab$BIC, digits),
                    CAIC = round(tab$CAIC, digits), ABIC = round(tab$ABIC, digits),
                    decision = ifelse(is.na(tab$decision), "", tab$decision))
  print(out, row.names = FALSE)
  for (ov in names(x$decisions))
    cat("\nExplanatory similarity with ", ov, " outcome variances: ", x$decisions[[ov]]$decision,
        " (", x$decisions[[ov]]$n_decrease, " of 3 criteria lower for equal means).\n", sep = "")
  nm <- names(x$fits)[1L]
  cat("\nOutcome means by profile (", nm, "):\n", sep = "")
  mm <- x$means[[nm]]
  mm$mean <- round(mm$mean, digits); mm$ci_lower <- round(mm$ci_lower, digits)
  mm$ci_upper <- round(mm$ci_upper, digits)
  print(mm[, c("group_level", "class", "variable", "mean", "ci_lower", "ci_upper", "proportion")],
        row.names = FALSE)
  print_estimation_notes(x$fits)
  invisible(x)
}
