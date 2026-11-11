#' Predictive similarity: predictors of profile membership across groups
#'
#' Tests whether the multinomial logistic relations between predictors and
#' profile membership are the same across groups (step 5 of Morin et al.,
#' 2016). The profile solution (indicator means and variances) is held fixed
#' at the estimates of a previously retained similarity model, as in the
#' start-value approach used with Mplus, and two models are compared: one with
#' group-specific predictor effects and one with effects constrained to be
#' equal across groups. Predictive similarity is supported when at least two
#' of CAIC, BIC and ABIC are lower for the constrained model.
#'
#' @param fit An `mglpa_fit` whose profile solution is held fixed, or an
#'   `mglpa_similarity` object (its most constrained supported model is used
#'   unless `source` names another one).
#' @param predictors Character vector of predictor names.
#' @param data Data frame; defaults to the data stored in `fit`.
#' @param source For an `mglpa_similarity` object, the name of the model whose
#'   solution is held fixed.
#' @param effects Which models to estimate: `"both"`, `"free"` or `"equal"`.
#' @param digits Rounding applied to the fixed parameter values (Mplus writes
#'   start values with five decimals); `NULL` keeps full precision.
#' @inheritParams fit_lpa
#' @return An object of class `mglpa_predictive` with elements `fits`
#'   (`free` and/or `equal`), `table` (fit indices and decision), `effects`
#'   (regression coefficients and odds ratios for all pairwise profile
#'   comparisons, by model) and `source` (the fixed profile solution).
#' @export
predictive_similarity <- function(fit, predictors, data = NULL, source = NULL,
                                  effects = c("both", "free", "equal"),
                                  starts = c(20, 5), stiterations = 10, maxiter = 500,
                                  tol = 1e-6, se = c("robust", "ml", "none"), digits = 5,
                                  seed = NULL, verbose = TRUE) {
  effects <- match.arg(effects)
  se <- match.arg(se)
  if (!is.null(seed)) set.seed(seed)
  src <- resolve_source(fit, source)
  if (is.null(data)) data <- src$data$data_frame
  if (is.null(data)) stop("`data` must be supplied.")
  src_fit <- src$fit
  m0 <- src_fit$model
  if (is.null(m0$group)) stop("Predictive similarity requires a model with a group variable.")
  fixed <- fixed_from_fit(src_fit, c("mean", "variance"), digits = digits)
  which_models <- if (effects == "both") c("free", "equal") else effects
  fits <- list()
  for (eff in which_models) {
    model <- mglpa_model(m0$k, m0$indicators, group = m0$group, variances = m0$variances,
                         equal_means = m0$equal_means, equal_variances = m0$equal_variances,
                         equal_proportions = m0$equal_proportions, predictors = predictors,
                         predictor_effects = eff, group_levels = m0$group_levels)
    f <- fit_lpa(data, model, starts = starts, stiterations = stiterations, maxiter = maxiter,
                 tol = tol, fixed = fixed, start_values = src_fit, se = se, verbose = FALSE)
    fits[[eff]] <- f
    if (verbose)
      message(sprintf("%-6s effects: LL = %.3f, npar = %d, BIC = %.3f, CAIC = %.3f, ABIC = %.3f",
                      eff, f$loglik, f$npar, f$fit[["BIC"]], f$fit[["CAIC"]], f$fit[["ABIC"]]))
  }
  rows <- lapply(names(fits), function(nm) fit_row(fits[[nm]], model = paste0(nm, "_effects")))
  tab <- do.call(rbind, rows)
  tab$group <- NULL
  tab$decision <- NA_character_
  tab$n_decrease <- NA_integer_
  decision <- NULL
  if (all(c("free", "equal") %in% names(fits))) {
    cmp <- compare_fit(fits$equal, fits$free)
    decision <- cmp
    i <- which(tab$model == "equal_effects")
    tab$decision[i] <- cmp$decision
    tab$n_decrease[i] <- cmp$n_decrease
  }
  eff_tables <- lapply(fits, predictor_effects)
  structure(list(fits = fits, table = tab, decision = decision, effects = eff_tables,
                 predictors = predictors, source = src$name, source_fit = src_fit,
                 data = data, call = match.call()),
            class = "mglpa_predictive")
}

resolve_source <- function(fit, source = NULL) {
  if (inherits(fit, "mglpa_similarity")) {
    nm <- if (is.null(source)) fit$retained else source
    if (!nm %in% names(fit$fits)) stop("Model `", nm, "` not found in the similarity tests.")
    return(list(fit = fit$fits[[nm]], name = nm, data = list(data_frame = fit$data)))
  }
  if (!inherits(fit, "mglpa_fit")) stop("`fit` must be an mglpa_fit or mglpa_similarity object.")
  list(fit = fit, name = similarity_level(fit$model), data = list(data_frame = fit$data_frame))
}

#' Multinomial regression effects of the predictors for all profile comparisons
#'
#' @param fit An `mglpa_fit` with predictors of profile membership.
#' @return A data frame with, for every group (or `"all"` when the effects
#'   are equal across groups), predictor and pair of profiles, the logit
#'   coefficient (profile vs. comparison profile), its standard error,
#'   z statistic, p-value, and the odds ratio with 95% confidence limits.
#' @export
predictor_effects <- function(fit) {
  stopifnot(inherits(fit, "mglpa_fit"))
  m <- fit$model
  if (!length(m$predictors)) stop("The model has no predictors.")
  K <- m$k
  G <- fit$D$G
  groups <- if (m$predictor_effects == "equal") list(list(g = NA, label = "all")) else
    lapply(seq_len(G), function(g) list(g = g, label = as.character(fit$data$group_levels[g])))
  L <- list()
  meta <- list()
  lab <- function(x, k, g) if (is.na(g)) sprintf("b_%s_c%d", x, k) else sprintf("b_%s_c%d_g%d", x, k, g)
  for (gr in groups) for (x in m$predictors) for (k1 in seq_len(K - 1L)) for (k2 in (k1 + 1L):K) {
    w <- c(1)
    names(w) <- lab(x, k1, gr$g)
    if (k2 < K) w <- c(w, stats::setNames(-1, lab(x, k2, gr$g)))
    id <- sprintf("%s|%s|%d|%d", gr$label, x, k1, k2)
    L[[id]] <- w
    meta[[id]] <- data.frame(group_level = gr$label, predictor = x, profile = k1,
                             reference = k2, comparison = sprintf("Profile %d vs. %d", k1, k2),
                             stringsAsFactors = FALSE)
  }
  ce <- contrast_estimates(fit, L, exponentiate = TRUE)
  out <- cbind(do.call(rbind, meta), ce[, c("estimate", "se", "z", "p", "exp_estimate",
                                            "exp_ci_lower", "exp_ci_upper")])
  names(out)[names(out) == "estimate"] <- "coef"
  names(out)[names(out) == "exp_estimate"] <- "OR"
  names(out)[names(out) == "exp_ci_lower"] <- "OR_ci_lower"
  names(out)[names(out) == "exp_ci_upper"] <- "OR_ci_upper"
  rownames(out) <- NULL
  out
}

#' @export
print.mglpa_predictive <- function(x, digits = 3, ...) {
  cat("<mglpa_predictive> predictors: ", paste(x$predictors, collapse = ", "),
      "; profile solution fixed at the ", x$source, " model\n", sep = "")
  tab <- x$table
  out <- data.frame(model = tab$model, LL = round(tab$LL, digits), npar = tab$npar,
                    AIC = round(tab$AIC, digits), BIC = round(tab$BIC, digits),
                    CAIC = round(tab$CAIC, digits), ABIC = round(tab$ABIC, digits),
                    decision = ifelse(is.na(tab$decision), "", tab$decision))
  print(out, row.names = FALSE)
  if (!is.null(x$decision))
    cat("\nPredictive similarity ", x$decision$decision, " (", x$decision$n_decrease,
        " of 3 information criteria lower for equal effects).\n", sep = "")
  for (nm in names(x$effects)) {
    cat("\n", nm, " effects (logit coefficients and odds ratios):\n", sep = "")
    e <- x$effects[[nm]]
    e$coef <- round(e$coef, digits); e$se <- round(e$se, digits); e$p <- round(e$p, 3)
    e$OR <- round(e$OR, digits)
    print(e[, c("group_level", "predictor", "comparison", "coef", "se", "p", "OR")], row.names = FALSE)
  }
  print_estimation_notes(x$fits)
  invisible(x)
}
