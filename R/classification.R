#' Classification quality of a profile solution
#'
#' Entropy and the classification tables that Mplus prints under
#' "Classification quality": the average posterior probabilities of profile
#' membership of the observations assigned to each profile (most likely
#' profile in rows, profiles in columns) and the probabilities of the most
#' likely profile given the latent profile (profiles in rows, most likely
#' profile in columns). For multi-group models the tables are computed within
#' each group.
#'
#' @param fit An `mglpa_fit` object.
#' @return An object of class `mglpa_classification`: a list with `entropy`
#'   (as reported in `fit`), `entropy_profiles`, `entropy_by_group`, `counts`
#'   (most likely profile counts by group), `average_posterior` and
#'   `classification_probabilities` (lists of matrices, one per group).
#' @examples
#' \donttest{
#' data(mglpa_wm)
#' ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
#' fit <- fit_lpa(mglpa_wm, mglpa_model(3, ind, group = "group", equal_means = TRUE),
#'                starts = c(20, 4), seed = 1)
#' classification_quality(fit)
#' }
#' @export
classification_quality <- function(fit) {
  stopifnot(inherits(fit, "mglpa_fit"))
  tau <- fit$posterior
  most <- fit$class
  K <- ncol(tau)
  D <- fit$D
  lev <- as.character(fit$data$group_levels)
  nm <- paste0("Profile ", seq_len(K))
  avg <- cp <- counts <- list()
  ent <- numeric(D$G)
  for (g in seq_len(D$G)) {
    rows <- D$grp_rows[[g]]
    tg <- tau[rows, , drop = FALSE]
    mg <- most[rows]
    A <- matrix(NA_real_, K, K, dimnames = list(nm, nm))
    P <- A
    n <- tabulate(mg, nbins = K)
    tot <- colSums(tg)
    for (c in seq_len(K)) {
      sel <- mg == c
      if (n[c] > 0) A[c, ] <- colMeans(tg[sel, , drop = FALSE])
      P[, c] <- colSums(tg[sel, , drop = FALSE]) / tot
    }
    avg[[lev[g]]] <- A
    cp[[lev[g]]] <- P
    counts[[lev[g]]] <- n
    ent[g] <- entropy_index(tg)
  }
  names(ent) <- lev
  structure(list(entropy = fit$entropy, entropy_profiles = fit$entropy_profiles,
                 entropy_by_group = ent, counts = counts, average_posterior = avg,
                 classification_probabilities = cp, k = K, group_levels = lev),
            class = "mglpa_classification")
}

#' @export
print.mglpa_classification <- function(x, digits = 3, ...) {
  multi <- length(x$group_levels) > 1L
  cat("Classification quality\n")
  if (is.na(x$entropy)) cat("  Entropy: not defined for one profile\n") else {
    cat(sprintf("  Entropy = %.*f", digits, x$entropy))
    if (multi) cat(sprintf(" (profiles x groups); profile-specific entropy = %.*f", digits,
                           x$entropy_profiles))
    cat("\n")
  }
  for (g in names(x$average_posterior)) {
    if (multi) cat(sprintf("\nGroup %s (entropy = %.*f)\n", g, digits, x$entropy_by_group[[g]]))
    cat("\nAverage posterior probabilities for most likely profile membership (row) by profile (column):\n")
    print(cbind(n = x$counts[[g]], round(x$average_posterior[[g]], digits)))
    cat("\nProbabilities of the most likely profile (column) given the profile (row):\n")
    print(round(x$classification_probabilities[[g]], digits))
  }
  invisible(x)
}

#' @rdname apa_tables
#' @export
apa_table.mglpa_classification <- function(x, digits = 3) {
  multi <- length(x$group_levels) > 1L
  rows <- list()
  for (g in names(x$average_posterior)) {
    A <- x$average_posterior[[g]]
    df <- data.frame(`Most likely profile` = seq_len(x$k), n = x$counts[[g]],
                     check.names = FALSE, stringsAsFactors = FALSE)
    for (j in seq_len(x$k)) df[[paste("Profile", j)]] <- fmt(A[, j], digits)
    if (multi) df <- cbind(Group = g, df, stringsAsFactors = FALSE)
    rows[[g]] <- df
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "note") <- paste0(
    "Average posterior probabilities of profile membership for the observations assigned to each ",
    "profile (most likely profile membership in rows, profiles in columns). ",
    if (is.na(x$entropy)) "" else sprintf("Entropy = %s.", fmt(x$entropy, digits)))
  out
}
