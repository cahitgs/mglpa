# Multinomial logistic part of the model: profile membership regressed on
# known-group dummies and (optionally) predictors. Coefficients are stored as a
# d x (K-1) matrix `th` whose columns correspond to the non-reference profiles;
# the last profile is the reference (as in Mplus).

logit_theta_from_par <- function(alpha, beta, info) {
  K <- info$K; G <- info$G; q <- info$q
  th <- matrix(0, ncol(info$Xd), K - 1L)
  for (g in seq_len(G)) {
    r <- info$int_row(g)
    th[r, ] <- alpha[, g]
    if (q > 0L) for (x in seq_len(q)) th[info$slp_row(x, g), ] <- beta[x, , g]
  }
  th
}

logit_par_from_theta <- function(th, info) {
  K <- info$K; G <- info$G; q <- info$q
  alpha <- array(0, c(K - 1L, G))
  beta <- if (q > 0L) array(0, c(q, K - 1L, G)) else NULL
  for (g in seq_len(G)) {
    alpha[, g] <- th[info$int_row(g), ]
    if (q > 0L) for (x in seq_len(q)) beta[x, , g] <- th[info$slp_row(x, g), ]
  }
  list(alpha = alpha, beta = beta)
}

softmax_ref <- function(eta) {
  # eta: N x (K-1); reference category has linear predictor 0
  N <- nrow(eta)
  full <- cbind(eta, 0)
  mx <- full[cbind(seq_len(N), max.col(full, ties.method = "first"))]
  lse <- mx + log(rowSums(exp(full - mx)))
  full - lse
}

# Weighted multinomial logistic regression with soft class labels `tau`,
# maximised by Newton-Raphson with step halving. Entries of `th` for which
# `mask` is FALSE are held fixed.
logit_newton <- function(tau, Xd, th, mask, maxit = 50L, tol = 1e-8) {
  d <- ncol(Xd); Km1 <- ncol(th); K <- Km1 + 1L
  f <- which(mask)
  if (!length(f)) return(th)
  Qfun <- function(th) sum(tau * softmax_ref(Xd %*% th))
  q0 <- Qfun(th)
  for (it in seq_len(maxit)) {
    logpi <- softmax_ref(Xd %*% th)
    pi <- exp(logpi)
    R <- tau[, -K, drop = FALSE] - pi[, -K, drop = FALSE]
    grad <- crossprod(Xd, R)
    gf <- grad[f]
    if (max(abs(gf)) < tol) break
    H <- matrix(0, d * Km1, d * Km1)
    for (k in seq_len(Km1)) for (m in k:Km1) {
      w <- if (k == m) pi[, k] * (1 - pi[, k]) else -pi[, k] * pi[, m]
      Hkm <- -crossprod(Xd * w, Xd)
      rk <- (k - 1L) * d + seq_len(d)
      cm <- (m - 1L) * d + seq_len(d)
      H[rk, cm] <- Hkm
      if (m != k) H[cm, rk] <- t(Hkm)
    }
    Hf <- H[f, f, drop = FALSE]
    step <- tryCatch(-solve(Hf, gf), error = function(e)
      -solve(Hf - diag(1e-6, length(f)), gf))
    s <- 1
    repeat {
      th_new <- th
      th_new[f] <- th[f] + s * step
      q1 <- Qfun(th_new)
      if (is.finite(q1) && q1 >= q0 - 1e-10) break
      s <- s / 2
      if (s < 1e-6) { th_new <- th; q1 <- q0; break }
    }
    if (q1 - q0 < 1e-12 && s < 1) { th <- th_new; break }
    th <- th_new
    q0 <- q1
  }
  th
}
