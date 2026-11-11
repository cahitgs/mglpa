#' @importFrom stats coef vcov logLik nobs
NULL

fit_indices <- function(ll, npar, N) {
  c(LL = ll, npar = npar, AIC = -2 * ll + 2 * npar, BIC = -2 * ll + npar * log(N),
    CAIC = -2 * ll + npar * (log(N) + 1), ABIC = -2 * ll + npar * log((N + 2) / 24))
}

# Entropy-based classification accuracy. With `n_patterns` the denominator
# uses the number of latent class patterns (profiles x known groups), which is
# what Mplus reports for KNOWNCLASS models; the posterior uncertainty in the
# numerator is identical because group membership is known.
entropy_index <- function(tau, n_patterns = ncol(tau)) {
  K <- ncol(tau)
  if (K < 2L) return(NA_real_)
  tl <- tau * log(tau)
  tl[!is.finite(tl)] <- 0
  1 - sum(-tl) / (nrow(tau) * log(n_patterns))
}

class_summary <- function(tau, es, D, K) {
  N <- nrow(tau)
  most <- max.col(tau, ties.method = "first")
  est_counts <- colSums(tau)
  ml_counts <- tabulate(most, nbins = K)
  pi <- exp(es$logpi)
  by_group <- function(f) matrix(unlist(lapply(seq_len(D$G), f)), nrow = D$G, ncol = K, byrow = TRUE)
  by_group_model <- by_group(function(g) colMeans(pi[D$grp_rows[[g]], , drop = FALSE]))
  by_group_posterior <- by_group(function(g) colMeans(tau[D$grp_rows[[g]], , drop = FALSE]))
  by_group_most <- by_group(function(g) tabulate(most[D$grp_rows[[g]]], nbins = K) /
                              length(D$grp_rows[[g]]))
  dimnames(by_group_model) <- dimnames(by_group_posterior) <- dimnames(by_group_most) <-
    list(as.character(D$group_levels), paste0("c", seq_len(K)))
  classif <- matrix(NA_real_, K, K, dimnames = list(paste0("c", seq_len(K)), paste0("c", seq_len(K))))
  for (k in seq_len(K)) if (ml_counts[k] > 0) classif[k, ] <- colMeans(tau[most == k, , drop = FALSE])
  list(
    estimated = data.frame(class = seq_len(K), count = est_counts, proportion = est_counts / N),
    most_likely = data.frame(class = seq_len(K), count = ml_counts, proportion = ml_counts / N),
    by_group_model = by_group_model,
    by_group_posterior = by_group_posterior,
    by_group_most_likely = by_group_most,
    classification = classif,
    most_likely_class = most,
    entropy = entropy_index(tau),
    entropy_patterns = entropy_index(tau, n_patterns = K * D$G)
  )
}

# ------------------------------------------------------------------------------
# Methods ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' @export
print.mglpa_fit <- function(x, digits = 3, ...) {
  cat("<mglpa_fit> ", model_description(x$model), "\n", sep = "")
  cat("  N = ", x$data$N, sep = "")
  if (length(x$data$n_g) > 1L)
    cat(" (", paste(sprintf("%s: %d", x$data$group_levels, x$data$n_g), collapse = ", "), ")", sep = "")
  cat("\n")
  cat(sprintf("  Log-likelihood = %.*f, free parameters = %d%s\n", digits, x$loglik, x$npar,
              if (x$mplus_convention && length(x$data$n_g) > 1L) " (Mplus KNOWNCLASS convention)" else ""))
  f <- x$fit
  cat(sprintf("  AIC = %.*f, BIC = %.*f, CAIC = %.*f, ABIC = %.*f\n", digits, f[["AIC"]],
              digits, f[["BIC"]], digits, f[["CAIC"]], digits, f[["ABIC"]]))
  if (!is.na(x$entropy)) cat(sprintf("  Entropy = %.*f\n", digits, x$entropy))
  st <- x$starts
  cat(sprintf("  Best log-likelihood replicated %d time(s) in %d final-stage optimisation(s)%s\n",
              st$best_ll_replications, length(st$final_ll),
              if (!x$converged) " [NOT CONVERGED]" else ""))
  pr <- x$class_counts$by_group_model
  cat("  Profile proportions (model-based):\n")
  print(round(pr, digits))
  print_estimation_notes(list(x))
  invisible(x)
}

# Mplus-style estimation warnings of a fit, built from its convergence status
# and the replication of its best log-likelihood over the final-stage
# optimisations. Used by fit_lpa_core() and, after relabelling the profiles,
# by reorder_profiles(), which keeps the random-start record of the original
# fit.
estimation_warnings <- function(converged, iterations, tol, K, n_final, reps, n_failed = 0L,
                                singular = FALSE) {
  warn <- character()
  if (!isTRUE(converged))
    warn <- c(warn, sprintf(paste(
      "THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY. The EM algorithm did not meet the",
      "convergence criterion (log-likelihood change < %g) within %d iterations. The results may",
      "not be trustworthy; increase `maxiter` or the number of random starts."), tol, iterations))
  # Replication of the best log-likelihood (only meaningful with several
  # final-stage optimisations and more than one profile)
  if (K > 1L && n_final >= 2L && reps == 1L)
    warn <- c(warn, paste(
      "THE BEST LOG-LIKELIHOOD VALUE WAS NOT REPLICATED. The solution may not be trustworthy",
      "due to local maxima. Increase the number of random starts."))
  else if (K > 1L && n_final >= 3L && reps < 3L)
    warn <- c(warn, sprintf(paste(
      "The best log-likelihood value was replicated only %d times (at least 3 recommended).",
      "Consider increasing the number of random starts."), reps))
  if (isTRUE(n_failed > 0L))
    warn <- c(warn, sprintf("%d of %d final-stage optimisations produced a non-finite log-likelihood.",
                            n_failed, n_final))
  if (isTRUE(singular))
    warn <- c(warn, paste(
      "THE STANDARD ERRORS MAY NOT BE TRUSTWORTHY: the information matrix is singular or",
      "ill-conditioned and a pseudo-inverse was used. This can be due to an empty or very small",
      "profile, a variance at its lower bound, or a non-identified model."))
  warn
}

# Print the estimation warnings of one or more fits (named list), if any.
print_estimation_notes <- function(fits) {
  has <- vapply(fits, function(f) length(f$warnings) > 0L, logical(1))
  if (!any(has)) return(invisible(NULL))
  cat("\nEstimation warnings:\n")
  for (nm in names(fits)[has]) for (w in fits[[nm]]$warnings)
    cat("  ", if (!is.null(nm) && nzchar(nm)) paste0(nm, ": ") else "", w, "\n", sep = "")
  invisible(NULL)
}

#' Parameter estimates of a fitted profile model
#'
#' @param object An `mglpa_fit` object.
#' @param type Optional parameter types to keep.
#' @param ... Unused.
#' @return A data frame with one row per model parameter (free or fixed).
#' @export
parameters <- function(object, type = NULL, ...) {
  stopifnot(inherits(object, "mglpa_fit"))
  pt <- object$params
  if (!is.null(type)) pt <- pt[pt$type %in% type, ]
  pt$group_level <- object$data$group_levels[pt$group]
  cols <- c("type", "variable", "class", "group", "group_level", "label", "free", "est", "se",
            "z", "p", "ci_lower", "ci_upper")
  out <- pt[, cols]
  rownames(out) <- NULL
  out
}

#' @export
summary.mglpa_fit <- function(object, ...) {
  structure(list(fit = object, parameters = parameters(object)), class = "summary.mglpa_fit")
}

#' @export
print.summary.mglpa_fit <- function(x, digits = 3, ...) {
  print(x$fit, digits = digits)
  cat("\nParameters (", x$fit$se_type, " standard errors):\n", sep = "")
  pt <- x$parameters
  pt <- pt[!duplicated(pt$label), ]
  pt$est <- round(pt$est, digits); pt$se <- round(pt$se, digits)
  pt$z <- round(pt$z, 2); pt$p <- round(pt$p, 3)
  print(pt[, c("type", "variable", "class", "group_level", "label", "free", "est", "se", "z", "p")],
        row.names = FALSE)
  invisible(x)
}

#' @export
coef.mglpa_fit <- function(object, ...) object$coefficients

#' @export
vcov.mglpa_fit <- function(object, ...) object$vcov

#' @export
logLik.mglpa_fit <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$data$N, class = "logLik")
}

#' @export
nobs.mglpa_fit <- function(object, ...) object$data$N

#' Posterior profile probabilities and most likely profile
#'
#' @param fit An `mglpa_fit` object.
#' @return A data frame with the posterior probabilities of each profile and
#'   the most likely profile for each observation used in the analysis.
#' @export
posterior_probabilities <- function(fit) {
  stopifnot(inherits(fit, "mglpa_fit"))
  out <- as.data.frame(fit$posterior)
  out$most_likely <- fit$class
  out
}
