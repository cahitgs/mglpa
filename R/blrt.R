#' Simulate data from a fitted single-group profile model
#'
#' Draws profile membership from the estimated profile proportions and the
#' indicators (and outcomes) from the profile-specific normal distributions.
#' Predictors of profile membership are not supported.
#'
#' @param fit An `mglpa_fit` estimated in a single group without predictors.
#' @param n Number of observations; defaults to the sample size of the fit.
#' @return A data frame with the indicator (and outcome) columns and a
#'   `profile` column with the generating profile.
#' @export
simulate_lpa <- function(fit, n = fit$data$N) {
  stopifnot(inherits(fit, "mglpa_fit"))
  if (fit$D$G > 1L || length(fit$model$predictors))
    stop("simulate_lpa() supports single-group models without predictors.")
  K <- fit$model$k
  p <- fit$D$p
  pt <- fit$params
  mu <- matrix(pt$est[pt$type == "mean"], nrow = p)
  s2 <- matrix(pt$est[pt$type == "variance"], nrow = p)
  pr <- fit$class_counts$by_group_model[1, ]
  cls <- sample.int(K, n, replace = TRUE, prob = pr)
  Y <- t(mu[, cls, drop = FALSE]) + matrix(stats::rnorm(n * p), n, p) * sqrt(t(s2[, cls, drop = FALSE]))
  out <- as.data.frame(Y)
  names(out) <- fit$D$vars
  out$profile <- cls
  out
}

#' Bootstrap likelihood ratio test for k versus k - 1 profiles
#'
#' Parametric bootstrap analogue of the Mplus `TECH14` output. Data sets are
#' generated from the (k - 1)-profile solution, the (k - 1)- and k-profile
#' models are estimated on each of them, and the bootstrap distribution of
#' the likelihood-ratio statistic 2(LL_k - LL_(k-1)) gives the p-value of
#' the observed statistic: the proportion of bootstrap statistics at least as
#' large as the observed one, as in Mplus. As in Mplus, the (k - 1)-profile
#' model is started from the parent estimates and the k-profile model uses
#' random starts.
#'
#' @details
#' For every draw the convergence of both models and the number of
#' final-stage optimisations of the k-profile model that replicated its best
#' log-likelihood are recorded (element `draws`). Two kinds of draws can
#' distort the p-value:
#' * the statistic of a draw is understated when the k-profile model stopped
#'   at a local maximum (best log-likelihood not replicated, estimation not
#'   converged, or a negative statistic, which is impossible at the global
#'   maximum of nested models). Such draws below the observed statistic could
#'   have exceeded it; like Mplus, `blrt()` reports them in a warning and
#'   gives the largest p-value they could produce (`p_upper`);
#' * the statistic of a draw is overstated when the (k - 1)-profile model did
#'   not converge. Such draws above the observed statistic might not exceed
#'   it; `p_lower` is the smallest p-value they could produce.
#'
#' Non-converged draws are included in the p-value (as in Mplus) and listed
#' in the warnings; draws in which no finite log-likelihood was obtained are
#' excluded and counted in `B - successful`. When `p`, `p_lower` and
#' `p_upper` lead to the same conclusion the result is robust to these
#' problems; otherwise increase the random starts of the k-profile model
#' (`starts`) or the iteration limit (`maxiter`).
#'
#' @param fit_k,fit_km1 `mglpa_fit` objects for k and k - 1 profiles
#'   estimated on the same single-group data with the same variance
#'   specification.
#' @param B Number of bootstrap draws.
#' @param starts Random starts for the k-profile model in each draw (initial
#'   stage, final stage); at least two final-stage optimisations are needed to
#'   judge the replication of the best log-likelihood.
#' @param seed Optional seed (for `cores > 1` the seed initialises the
#'   parallel random number streams of the worker processes).
#' @param verbose Logical; report progress.
#' @param cores Number of CPU cores over which the bootstrap draws are
#'   distributed (needs the installed package, see [fit_lpa()]).
#' @param maxiter,tol Iteration limit and convergence criterion of the EM
#'   algorithm in the bootstrap draws; by default those of `fit_k`.
#' @return A list of class `mglpa_blrt` with the observed `statistic`, the
#'   `bootstrap` statistics, a data frame `draws` with the per-draw statistic,
#'   convergence and replication record, the p-value `p` and its bounds
#'   `p_lower` and `p_upper`, the number of draws requested `B`, the number
#'   of `successful` draws, the numbers of draws flagged as understated
#'   (`n_flagged`), overstated (`n_flagged_high`) and non-converged
#'   (`n_not_converged`), and the `warnings`.
#' @export
blrt <- function(fit_k, fit_km1, B = 100, starts = c(20, 4), seed = NULL, verbose = TRUE,
                 cores = 1L, maxiter = NULL, tol = NULL) {
  stopifnot(inherits(fit_k, "mglpa_fit"), inherits(fit_km1, "mglpa_fit"))
  if (fit_k$model$k != fit_km1$model$k + 1L) stop("`fit_k` must have one more profile than `fit_km1`.")
  if (fit_k$data$N != fit_km1$data$N) stop("The two fits must be estimated on the same data.")
  # Evaluate the arguments here: the function that analyses a draw is sent to
  # the worker processes, where an argument given as a variable of the calling
  # environment could not be evaluated any more.
  starts <- as.integer(starts)
  if (length(starts) != 2L || anyNA(starts) || any(starts < 0L))
    stop("`starts` must be two non-negative integers.")
  verbose <- isTRUE(verbose)
  cores <- max(1L, as.integer(cores))
  B <- as.integer(B)
  if (is.null(maxiter)) maxiter <- if (is.null(fit_k$starts$maxiter)) 500 else fit_k$starts$maxiter
  if (is.null(tol)) tol <- if (is.null(fit_k$starts$tol)) 1e-6 else fit_k$starts$tol
  if (!is.null(seed)) set.seed(seed)
  stat_obs <- 2 * (fit_k$loglik_mixture - fit_km1$loglik_mixture)
  m_km1 <- fit_km1$model
  m_k <- fit_k$model
  failed <- c(statistic = NA_real_, converged_k = NA_real_, converged_km1 = NA_real_,
              replications = NA_real_, n_final = NA_real_)
  one_draw <- function(b) {
    d <- simulate_lpa(fit_km1)
    tryCatch(suppressWarnings({
      f0 <- fit_lpa(d, m_km1, starts = c(0, 0), start_values = fit_km1, se = "none",
                    maxiter = maxiter, tol = tol, verbose = FALSE)
      f1 <- fit_lpa(d, m_k, starts = starts, se = "none", maxiter = maxiter, tol = tol,
                    verbose = FALSE)
      c(statistic = 2 * (f1$loglik_mixture - f0$loglik_mixture),
        converged_k = as.numeric(isTRUE(f1$converged)),
        converged_km1 = as.numeric(isTRUE(f0$converged)),
        replications = as.numeric(f1$starts$best_ll_replications),
        n_final = as.numeric(length(f1$starts$final_ll)))
    }), error = function(e) failed)
  }
  cl <- NULL
  if (cores > 1L && B > 1L) {
    cl <- start_cluster(cores)
    if (!is.null(cl)) {
      on.exit(parallel::stopCluster(cl), add = TRUE)
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
    }
  }
  chunk <- blrt_chunk(if (is.null(cl)) 1L else cores)
  rows <- list()
  for (idx in split(seq_len(B), ceiling(seq_len(B) / chunk))) {
    res <- if (is.null(cl)) lapply(idx, one_draw) else parallel::parLapply(cl, idx, one_draw)
    rows <- c(rows, res)
    if (verbose)
      message(sprintf("BLRT %d vs %d profiles: draw %d/%d", m_k$k, m_km1$k, length(rows), B))
  }
  draws <- as.data.frame(do.call(rbind, rows))
  draws$converged_k <- draws$converged_k == 1
  draws$converged_km1 <- draws$converged_km1 == 1
  summ <- blrt_summary(draws, stat_obs)
  if (length(summ$warnings))
    warning(paste(c(sprintf("mglpa BLRT (%d vs %d profiles):", m_k$k, m_km1$k), summ$warnings),
                  collapse = "\n  "), call. = FALSE)
  structure(c(list(k = m_k$k, statistic = stat_obs, bootstrap = draws$statistic, draws = draws,
                   B = B, starts = starts), summ),
            class = "mglpa_blrt")
}

# Draws are processed (and progress is reported) in chunks of this size
blrt_chunk <- function(cores) if (cores > 1L) max(10L, as.integer(cores)) else 10L

# p-value and local-maxima / convergence diagnostics of the bootstrap draws
# (see blrt()); `draws` has the columns statistic, converged_k, converged_km1,
# replications and n_final.
blrt_summary <- function(draws, stat_obs) {
  ok <- is.finite(draws$statistic)
  n_ok <- sum(ok)
  exceed <- ok & draws$statistic >= stat_obs
  n_exceed <- sum(exceed)
  p <- if (n_ok) n_exceed / n_ok else NA_real_
  # Understated statistic: the k-profile model stopped at a local maximum (not
  # replicated, not converged) or the statistic is negative, which is
  # impossible at the global maximum of nested models.
  understated <- ok & (!draws$converged_k | (draws$n_final >= 2 & draws$replications < 2) |
                         draws$statistic < 0)
  flagged_low <- understated & !exceed
  # Overstated statistic: the (k - 1)-profile model did not converge.
  flagged_high <- ok & !draws$converged_km1 & exceed
  n_low <- sum(flagged_low)
  n_high <- sum(flagged_high)
  p_upper <- if (n_ok) (n_exceed + n_low) / n_ok else NA_real_
  p_lower <- if (n_ok) (n_exceed - n_high) / n_ok else NA_real_
  nc_k <- sum(ok & !draws$converged_k)
  nc_km1 <- sum(ok & !draws$converged_km1)
  n_nc <- sum(ok & (!draws$converged_k | !draws$converged_km1))
  warn <- character()
  if (n_low > 0L)
    warn <- c(warn, sprintf(paste(
      "OF THE %d BOOTSTRAP DRAWS, %d HAD BOTH A SMALLER LRT VALUE THAN THE OBSERVED LRT VALUE AND",
      "NOT A REPLICATED (OR NOT CONVERGED) BEST LOG-LIKELIHOOD VALUE. THE P-VALUE MAY NOT BE",
      "TRUSTWORTHY DUE TO LOCAL MAXIMA (it could be as large as %.4f). Increase the number of",
      "random starts of the k-profile model (`starts`)."), n_ok, n_low, p_upper))
  if (n_high > 0L)
    warn <- c(warn, sprintf(paste(
      "%d of the %d draws with an LRT value at least as large as the observed one had a",
      "(k - 1)-profile model that did not converge, so their LRT value may be overstated;",
      "the p-value could be as small as %.4f. Increase `maxiter`."), n_high, n_exceed, p_lower))
  if (n_nc > 0L)
    warn <- c(warn, sprintf(paste(
      "%d of %d bootstrap draws did not meet the convergence criterion (k-profile model: %d,",
      "(k - 1)-profile model: %d); they are included in the p-value. Consider increasing `maxiter`."),
      n_nc, n_ok, nc_k, nc_km1))
  if (n_ok < nrow(draws))
    warn <- c(warn, sprintf("%d of %d bootstrap draws failed (no finite log-likelihood) and were excluded.",
                            nrow(draws) - n_ok, nrow(draws)))
  list(p = p, p_lower = p_lower, p_upper = p_upper, successful = n_ok, n_flagged = n_low,
       n_flagged_high = n_high, n_not_converged = n_nc, warnings = warn)
}

#' @export
print.mglpa_blrt <- function(x, ...) {
  cat(sprintf(paste("Bootstrap LRT for %d vs %d profiles: 2*(LL_k - LL_k-1) = %.3f,",
                    "p = %.4f (%d successful of %d draws)\n"),
              x$k, x$k - 1L, x$statistic, x$p, x$successful, x$B))
  if (isTRUE(x$n_flagged > 0L) || isTRUE(x$n_flagged_high > 0L))
    cat(sprintf("  Possible local maxima or non-convergence in %d draw(s): p could range from %.4f to %.4f\n",
                x$n_flagged + x$n_flagged_high, x$p_lower, x$p_upper))
  for (w in x$warnings) cat("  Warning: ", w, "\n", sep = "")
  invisible(x)
}

#' Add bootstrap likelihood ratio tests to a class enumeration
#'
#' Runs [blrt()] for every consecutive pair of solutions in each group and
#' adds a `BLRT_p` column to the enumeration table; the warnings of the tests
#' are appended to the `notes` column.
#'
#' @param enum An `mglpa_enumeration` object.
#' @inheritParams blrt
#' @return The updated `mglpa_enumeration` object, with the test objects in
#'   the element `blrt`.
#' @export
enumeration_blrt <- function(enum, B = 100, starts = c(20, 4), seed = NULL, verbose = TRUE,
                             cores = 1L, maxiter = NULL, tol = NULL) {
  stopifnot(inherits(enum, "mglpa_enumeration"))
  if (!is.null(seed)) set.seed(seed)
  enum$table$BLRT_p <- NA_real_
  enum$blrt <- list()
  for (g in names(enum$fits)) {
    ks <- as.integer(names(enum$fits[[g]]))
    for (k in ks[-1]) {
      if (!(k - 1L) %in% ks) next
      res <- blrt(enum$fits[[g]][[as.character(k)]], enum$fits[[g]][[as.character(k - 1L)]],
                  B = B, starts = starts, verbose = verbose, cores = cores,
                  maxiter = maxiter, tol = tol,
                  seed = if (is.null(seed)) NULL else seed + k)
      enum$blrt[[paste(g, k)]] <- res
      i <- enum$table$group == g & enum$table$k == k
      enum$table$BLRT_p[i] <- res$p
      if (length(res$warnings))
        enum$table$notes[i] <- trimws(paste(enum$table$notes[i],
                                            paste0("BLRT: ", res$warnings, collapse = " ")))
    }
  }
  enum$table <- order_enum_columns(enum$table)
  enum
}
