#' Specify a (multi-group) latent profile model
#'
#' Describes a latent profile analysis (LPA) model with continuous indicators
#' that is estimated in one sample or across known groups. The equality
#' constraints across groups correspond to the similarity tests of Morin,
#' Meyer, Creusier and Bietry (2016): the configural model imposes none, the
#' structural model adds equal profile means, the dispersion model adds equal
#' within-profile variances and the distributional model adds equal profile
#' proportions. The model can be extended with predictors of profile
#' membership (multinomial logistic regression; predictive similarity) and with
#' distal outcomes modelled as additional indicators (explanatory similarity).
#'
#' @param k Number of latent profiles.
#' @param indicators Character vector naming the profile indicators.
#' @param group Optional name of the known-group variable. When `NULL` the
#'   model is estimated in a single group.
#' @param variances `"class_invariant"` keeps the within-profile variances
#'   equal across profiles (the Mplus default); `"class_varying"` estimates a
#'   separate variance in each profile.
#' @param equal_means,equal_variances,equal_proportions Logical. Constrain the
#'   indicator means, the indicator variances or the profile proportions to be
#'   equal across groups.
#' @param predictors Optional character vector of predictors of profile
#'   membership.
#' @param predictor_effects `"free"` estimates group-specific multinomial
#'   regression coefficients; `"equal"` constrains them to equality across
#'   groups.
#' @param outcomes Optional character vector of distal outcomes.
#' @param outcome_means `"free"` estimates group- and profile-specific outcome
#'   means; `"equal"` constrains them to equality across groups.
#' @param outcome_variances `"equal"` estimates one residual variance per
#'   outcome; `"free"` estimates a profile- and group-specific variance.
#' @param group_levels Optional vector fixing the order of the group codes.
#'   Defaults to the sorted unique values found in the data.
#' @return An object of class `mglpa_model`.
#' @references Morin, A. J. S., Meyer, J. P., Creusier, J., & Bietry, F.
#'   (2016). Multiple-group analysis of similarity in latent profile solutions.
#'   *Organizational Research Methods, 19*(2), 231-254.
#' @examples
#' mglpa_model(3, c("verbal_wm", "visual_wm", "updating"), group = "group",
#'             equal_means = TRUE)
#' @export
mglpa_model <- function(k, indicators, group = NULL,
                        variances = c("class_invariant", "class_varying"),
                        equal_means = FALSE, equal_variances = FALSE,
                        equal_proportions = FALSE,
                        predictors = NULL,
                        predictor_effects = c("free", "equal"),
                        outcomes = NULL,
                        outcome_means = c("free", "equal"),
                        outcome_variances = c("equal", "free"),
                        group_levels = NULL) {
  variances <- match.arg(variances)
  predictor_effects <- match.arg(predictor_effects)
  outcome_means <- match.arg(outcome_means)
  outcome_variances <- match.arg(outcome_variances)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || k != round(k))
    stop("`k` must be a single positive integer.", call. = FALSE)
  if (!is.character(indicators) || !length(indicators))
    stop("`indicators` must be a non-empty character vector.", call. = FALSE)
  all_vars <- c(indicators, outcomes, predictors, group)
  if (anyDuplicated(all_vars))
    stop("A variable can appear only once across indicators, outcomes, predictors and group.",
         call. = FALSE)
  if (is.null(group) && (equal_means || equal_variances || equal_proportions))
    stop("Equality constraints across groups require a `group` variable.", call. = FALSE)
  if (!is.null(group) && (!is.character(group) || length(group) != 1L))
    stop("`group` must be a single variable name.", call. = FALSE)
  structure(list(
    k = as.integer(k), indicators = indicators, group = group,
    variances = variances, equal_means = isTRUE(equal_means),
    equal_variances = isTRUE(equal_variances),
    equal_proportions = isTRUE(equal_proportions),
    predictors = predictors, predictor_effects = predictor_effects,
    outcomes = outcomes, outcome_means = outcome_means,
    outcome_variances = outcome_variances, group_levels = group_levels
  ), class = "mglpa_model")
}

#' @export
print.mglpa_model <- function(x, ...) {
  cat("<mglpa_model> ", model_description(x), "\n", sep = "")
  cat("  indicators: ", paste(x$indicators, collapse = ", "), "\n", sep = "")
  if (!is.null(x$group)) cat("  group: ", x$group, "\n", sep = "")
  if (length(x$predictors))
    cat("  predictors: ", paste(x$predictors, collapse = ", "),
        " (", x$predictor_effects, " effects)\n", sep = "")
  if (length(x$outcomes))
    cat("  outcomes: ", paste(x$outcomes, collapse = ", "),
        " (", x$outcome_means, " means, ", x$outcome_variances, " variances)\n", sep = "")
  invisible(x)
}

model_description <- function(model) {
  parts <- sprintf("%d-profile LPA", model$k)
  parts <- c(parts, if (model$variances == "class_invariant") "class-invariant variances"
             else "class-varying variances")
  if (!is.null(model$group)) {
    cons <- c(if (model$equal_means) "means", if (model$equal_variances) "variances",
              if (model$equal_proportions) "proportions")
    parts <- c(parts, if (length(cons)) paste("equal across groups:", paste(cons, collapse = ", "))
               else "no cross-group constraints")
  }
  paste(parts, collapse = "; ")
}

similarity_level <- function(model) {
  if (is.null(model$group)) return("single group")
  if (model$equal_proportions && model$equal_variances && model$equal_means) return("distributional")
  if (model$equal_variances && model$equal_means) return("dispersion")
  if (model$equal_means) return("structural")
  "configural"
}

# ------------------------------------------------------------------------------
# Parameter table ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# One row per model parameter (free or fixed). Rows that share a label are a
# single free parameter (Mplus-style equality constraints). Rows are created in
# array order (variable fastest, then class, then group) so that the values of
# each type can be reshaped directly into p x K x G arrays.
build_param_table <- function(model, G) {
  K <- model$k
  ind <- model$indicators
  out <- model$outcomes
  pred <- model$predictors
  vars <- c(ind, out)
  p <- length(vars)
  q <- length(pred)
  n_ind <- length(ind)

  mean_label <- function(v, is_out, k, g) {
    equal <- if (is_out) model$outcome_means == "equal" else model$equal_means
    if (equal) sprintf("m_%s_c%d", v, k) else sprintf("m_%s_c%d_g%d", v, k, g)
  }
  var_label <- function(v, is_out, k, g) {
    if (is_out) {
      if (model$outcome_variances == "free") sprintf("v_%s_c%d_g%d", v, k, g)
      else sprintf("v_%s", v)
    } else {
      cv <- model$variances == "class_varying"
      eg <- model$equal_variances
      if (cv && !eg) sprintf("v_%s_c%d_g%d", v, k, g)
      else if (cv && eg) sprintf("v_%s_c%d", v, k)
      else if (!cv && !eg) sprintf("v_%s_g%d", v, g)
      else sprintf("v_%s", v)
    }
  }

  rows <- list()
  add <- function(type, variable, var_index, class, group, label, pos) {
    rows[[length(rows) + 1L]] <<- list(type = type, variable = variable,
                                       var_index = var_index, class = class,
                                       group = group, label = label, pos = pos)
  }
  for (g in seq_len(G)) for (k in seq_len(K)) for (j in seq_len(p)) {
    add("mean", vars[j], j, k, g, mean_label(vars[j], j > n_ind, k, g),
        (g - 1L) * p * K + (k - 1L) * p + j)
  }
  for (g in seq_len(G)) for (k in seq_len(K)) for (j in seq_len(p)) {
    add("variance", vars[j], j, k, g, var_label(vars[j], j > n_ind, k, g),
        (g - 1L) * p * K + (k - 1L) * p + j)
  }
  if (K > 1L) {
    for (g in seq_len(G)) for (k in seq_len(K - 1L)) {
      add("intercept", "", NA_integer_, k, g,
          if (model$equal_proportions) sprintf("a_c%d", k) else sprintf("a_c%d_g%d", k, g),
          (g - 1L) * (K - 1L) + k)
    }
    if (q > 0L) {
      for (g in seq_len(G)) for (k in seq_len(K - 1L)) for (x in seq_len(q)) {
        add("slope", pred[x], x, k, g,
            if (model$predictor_effects == "equal") sprintf("b_%s_c%d", pred[x], k)
            else sprintf("b_%s_c%d_g%d", pred[x], k, g),
            (g - 1L) * q * (K - 1L) + (k - 1L) * q + x)
      }
    }
  }
  pt <- data.frame(
    type = vapply(rows, `[[`, "", "type"),
    variable = vapply(rows, `[[`, "", "variable"),
    var_index = vapply(rows, `[[`, 1L, "var_index"),
    class = vapply(rows, `[[`, 1L, "class"),
    group = vapply(rows, `[[`, 1L, "group"),
    label = vapply(rows, `[[`, "", "label"),
    pos = vapply(rows, `[[`, 1L, "pos"),
    stringsAsFactors = FALSE
  )
  pt$free <- TRUE
  pt$value <- NA_real_
  pt
}

apply_fixed <- function(pt, fixed) {
  if (is.null(fixed) || !length(fixed)) return(pt)
  if (is.null(names(fixed)) || !is.numeric(fixed))
    stop("`fixed` must be a named numeric vector of parameter labels.", call. = FALSE)
  unknown <- setdiff(names(fixed), pt$label)
  if (length(unknown))
    stop("Unknown parameter labels in `fixed`: ", paste(unknown, collapse = ", "), call. = FALSE)
  hit <- pt$label %in% names(fixed)
  pt$free[hit] <- FALSE
  pt$value[hit] <- fixed[pt$label[hit]]
  pt
}

#' Extract parameter values from a fit to fix them in a later model
#'
#' Returns a named vector of estimates (named by parameter label) that can be
#' passed to the `fixed` argument of [fit_lpa()]. This is how the predictive
#' and explanatory similarity models hold the profile solution constant at the
#' values obtained in an earlier step, in the spirit of the start-value
#' approach used with Mplus.
#'
#' @param fit An `mglpa_fit` object.
#' @param types Parameter types to extract: any of `"mean"`, `"variance"`,
#'   `"intercept"` (profile-membership logits) and `"slope"`.
#' @param digits Optional rounding of the extracted values (Mplus writes start
#'   values with five decimals). `NULL` keeps full precision.
#' @return A named numeric vector.
#' @export
fixed_from_fit <- function(fit, types = c("mean", "variance"), digits = NULL) {
  stopifnot(inherits(fit, "mglpa_fit"))
  pt <- fit$params
  pt <- pt[pt$type %in% types, ]
  vals <- pt$est
  names(vals) <- pt$label
  vals <- vals[!duplicated(names(vals))]
  if (!is.null(digits)) vals <- round(vals, digits)
  vals
}
