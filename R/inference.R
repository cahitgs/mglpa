# Analytic per-observation score vectors of the observed-data log-likelihood
# with respect to the free parameters, posterior-weighted complete-data scores
# (Louis, 1982). Returns an N x d matrix.
score_matrix <- function(par, es, D, info) {
  pt <- info$pt
  N <- D$N; K <- info$K
  S <- matrix(0, N, info$d)
  tau <- es$tau
  pi <- if (K > 1L) exp(es$logpi) else NULL
  for (r in which(pt$free)) {
    type <- pt$type[r]; k <- pt$class[r]; g <- pt$group[r]; col <- pt$index[r]
    rows <- D$grp_rows[[g]]
    if (type == "mean" || type == "variance") {
      j <- pt$var_index[r]
      y <- D$Y[rows, j]
      o <- D$O[rows, j]
      mu <- par$mu[j, k, g]; s2 <- par$sig2[j, k, g]
      contrib <- if (type == "mean") tau[rows, k] * (y - mu) / s2 else
        tau[rows, k] * (-0.5 / s2 + 0.5 * (y - mu)^2 / (s2 * s2))
      contrib[!o] <- 0
    } else if (type == "intercept") {
      contrib <- tau[rows, k] - pi[rows, k]
    } else {
      x <- pt$var_index[r]
      contrib <- (tau[rows, k] - pi[rows, k]) * D$grp[[g]]$X[, x]
    }
    S[rows, col] <- S[rows, col] + contrib
  }
  S
}

loglik_gradient <- function(theta, D, info) {
  par <- unpack_theta(theta, info)
  es <- estep(par, D, info)
  colSums(score_matrix(par, es, D, info))
}

# Observed information by central differences of the analytic gradient,
# its inverse (ML covariance matrix) and the per-observation scores at the MLE.
observed_information <- function(theta, pt, D, info, model) {
  d <- length(theta)
  par <- unpack_theta(theta, info)
  es <- estep(par, D, info)
  S <- score_matrix(par, es, D, info)
  H <- matrix(0, d, d)
  for (m in seq_len(d)) {
    h <- 1e-5 * max(1, abs(theta[m]))
    tp <- theta; tp[m] <- tp[m] + h
    tm <- theta; tm[m] <- tm[m] - h
    H[, m] <- (loglik_gradient(tp, D, info) - loglik_gradient(tm, D, info)) / (2 * h)
  }
  H <- (H + t(H)) / 2
  I <- -H
  ev <- tryCatch(eigen(I, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA)
  cond <- if (all(is.finite(ev)) && max(ev) > 0) min(ev) / max(ev) else NA_real_
  V <- safe_inverse(I)
  singular <- isTRUE(attr(V, "singular")) || (is.finite(cond) && cond < 1e-12) ||
    any(!is.finite(ev)) || any(ev <= 0)
  attr(V, "singular") <- NULL
  list(information = I, vcov = V, scores = S, condition = cond, gradient = colSums(S),
       singular = singular)
}

# Inverse of the information matrix, or a pseudo-inverse (flagged with the
# attribute "singular") when it cannot be inverted.
safe_inverse <- function(M) {
  out <- tryCatch(solve(M), error = function(e) NULL)
  singular <- FALSE
  if (is.null(out) || any(!is.finite(out))) {
    singular <- TRUE
    s <- svd(M)
    pos <- s$d > max(s$d) * 1e-10
    out <- s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
  }
  attr(out, "singular") <- singular
  out
}

#' Linear combinations of parameters with delta-method standard errors
#'
#' @param fit An `mglpa_fit` object.
#' @param L A matrix with one row per contrast and one column per free
#'   parameter (named by label), or a named list of named numeric vectors
#'   giving the weights of the labels involved in each contrast.
#' @param exponentiate Logical; also report `exp()` of the estimates and
#'   confidence limits (odds ratios).
#' @return A data frame with estimates, standard errors, z statistics,
#'   p-values and 95% confidence limits.
#' @export
contrast_estimates <- function(fit, L, exponentiate = FALSE) {
  stopifnot(inherits(fit, "mglpa_fit"))
  if (is.null(fit$vcov)) stop("The fit has no covariance matrix (se = 'none').")
  theta <- fit$coefficients
  if (is.list(L)) {
    Lm <- matrix(0, length(L), length(theta), dimnames = list(names(L), names(theta)))
    for (i in seq_along(L)) {
      w <- L[[i]]
      unknown <- setdiff(names(w), names(theta))
      if (length(unknown)) stop("Unknown parameter labels: ", paste(unknown, collapse = ", "))
      Lm[i, names(w)] <- w
    }
    L <- Lm
  }
  est <- as.vector(L %*% theta)
  V <- L %*% fit$vcov %*% t(L)
  se <- sqrt(pmax(diag(V), 0))
  z <- est / se
  out <- data.frame(contrast = rownames(L), estimate = est, se = se, z = z,
                    p = 2 * stats::pnorm(-abs(z)),
                    ci_lower = est - stats::qnorm(0.975) * se,
                    ci_upper = est + stats::qnorm(0.975) * se,
                    stringsAsFactors = FALSE)
  if (exponentiate) {
    out$exp_estimate <- exp(out$estimate)
    out$exp_ci_lower <- exp(out$ci_lower)
    out$exp_ci_upper <- exp(out$ci_upper)
  }
  rownames(out) <- NULL
  out
}
