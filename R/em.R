LOG2PI <- log(2 * pi)

#' Fit a latent profile model by maximum likelihood
#'
#' Estimates the model described by [mglpa_model()] with an
#' expectation-maximisation (EM) algorithm. Profile indicators and distal
#' outcomes are conditionally independent normal variables given the profile
#' (a diagonal Gaussian mixture), profile membership follows a multinomial
#' logistic model with optional predictors and known groups, and cross-group
#' equality constraints are handled by pooling the sufficient statistics of
#' the parameters that share a label. Missing indicator values are handled by
#' full-information maximum likelihood (each observation contributes the
#' density of its observed indicators). The procedure mirrors the random-start
#' strategy of Mplus: `starts[1]` random starts are run for `stiterations`
#' iterations, the best `starts[2]` are iterated to convergence and the
#' highest log-likelihood is retained.
#'
#' @param data A data frame.
#' @param model An `mglpa_model` object.
#' @param starts Integer vector of length two: number of initial-stage random
#'   starts and number of final-stage optimisations (like `STARTS` in Mplus).
#' @param stiterations Number of EM iterations in the initial stage.
#' @param maxiter Maximum number of EM iterations in the final stage.
#' @param tol Convergence criterion on the absolute change in log-likelihood.
#' @param fixed Optional named numeric vector of parameter labels and the
#'   values at which they are fixed; see [fixed_from_fit()].
#' @param start_values Optional named numeric vector (parameter labels) or
#'   `mglpa_fit` whose estimates are used as an additional, unperturbed start.
#' @param se Type of standard errors: `"robust"` (sandwich estimator with the
#'   observed information matrix, equivalent to `ESTIMATOR = MLR`), `"ml"`
#'   (inverse observed information) or `"none"`.
#' @param mplus_convention Logical. When `TRUE` (the default) the
#'   log-likelihood includes the contribution of the known-group proportions
#'   and these proportions are counted as free parameters, which reproduces
#'   the log-likelihood and information criteria reported by Mplus for
#'   `KNOWNCLASS` models.
#' @param seed Optional random seed.
#' @param verbose Logical; print progress.
#' @param engine `"auto"` runs the EM iterations in compiled C++ code
#'   whenever the model allows it (no predictors of profile membership and
#'   free membership logits) and in R otherwise; `"r"` forces the pure-R
#'   implementation, which gives identical results and is kept as a
#'   reference.
#' @param cores Number of CPU cores used to distribute the random starts
#'   (a PSOCK cluster from the `parallel` package). Requires the installed
#'   package; falls back to sequential estimation with a warning otherwise.
#'   Results depend on the random starts, not on the number of cores, but the
#'   assignment of starts to workers makes exact reproducibility require the
#'   same `seed` and `cores`.
#' @return An object of class `mglpa_fit`.
#' @examples
#' \donttest{
#' data(mglpa_wm)
#' m <- mglpa_model(3, c("verbal_wm", "visual_wm", "updating", "inhibition",
#'                       "proc_speed"))
#' fit <- fit_lpa(mglpa_wm[mglpa_wm$group == 1, ], m, starts = c(20, 4), seed = 1)
#' fit
#' }
#' @export
fit_lpa <- function(data, model, starts = c(100, 10), stiterations = 10,
                    maxiter = 500, tol = 1e-6, fixed = NULL,
                    start_values = NULL, se = c("robust", "ml", "none"),
                    mplus_convention = TRUE, seed = NULL, verbose = FALSE,
                    engine = c("auto", "r"), cores = 1L) {
  se <- match.arg(se)
  engine <- match.arg(engine)
  stopifnot(inherits(model, "mglpa_model"))
  D <- prepare_data(data, model)
  fit_lpa_core(D, model, starts = starts, stiterations = stiterations, maxiter = maxiter,
               tol = tol, fixed = fixed, start_values = start_values, se = se,
               mplus_convention = mplus_convention, seed = seed, verbose = verbose,
               data_frame = data, engine = engine, cores = cores)
}

# Core estimation working on prepared data (see prepare_data()). Used by
# fit_lpa() and by functions that refit an existing model, e.g. reorder_profiles().
fit_lpa_core <- function(D, model, starts = c(100, 10), stiterations = 10, maxiter = 500,
                         tol = 1e-6, fixed = NULL, start_values = NULL, se = "robust",
                         mplus_convention = TRUE, seed = NULL, verbose = FALSE,
                         data_frame = NULL, engine = "auto", cores = 1L) {
  if (length(starts) != 2L || any(starts < 0)) stop("`starts` must be two non-negative integers.")
  if (!is.null(seed)) set.seed(seed)
  t0 <- Sys.time()

  pt <- build_param_table(model, D$G)
  pt <- apply_fixed(pt, fixed)
  info <- make_info(pt, D, model, engine = engine)
  pt <- info$pt
  if (info$d == 0L) stop("The model has no free parameters.")

  user_vals <- NULL
  if (!is.null(start_values)) {
    sv <- if (inherits(start_values, "mglpa_fit")) stats::coef(start_values) else start_values
    user_vals <- random_start(pt, D, info, model)
    hit <- pt$free & pt$label %in% names(sv)
    user_vals[hit] <- sv[pt$label[hit]]
  }

  res <- run_random_starts(pt, D, info, model, starts, stiterations, maxiter, tol,
                           user_vals, verbose, cores = cores)
  best <- res$best
  par <- best$par
  theta <- pack_theta(par, info)
  es <- best$es

  # Parameter table with estimates --------------------------------------------
  vals <- pt$value
  vals[pt$free] <- theta[pt$index[pt$free]]
  pt$est <- vals
  pt$se <- NA_real_

  vc <- NULL; vc_ml <- NULL; cond <- NA_real_; S <- NULL; singular <- FALSE
  if (se != "none") {
    inf <- observed_information(theta, pt, D, info, model)
    S <- inf$scores
    vc_ml <- inf$vcov
    cond <- inf$condition
    singular <- inf$singular
    vc <- if (se == "robust") vc_ml %*% crossprod(S) %*% vc_ml else vc_ml
    dimnames(vc) <- list(info$free_labels, info$free_labels)
    sev <- sqrt(pmax(diag(vc), 0))
    pt$se[pt$free] <- sev[pt$index[pt$free]]
  }
  pt$z <- pt$est / pt$se
  pt$p <- 2 * stats::pnorm(-abs(pt$z))
  pt$ci_lower <- pt$est - stats::qnorm(0.975) * pt$se
  pt$ci_upper <- pt$est + stats::qnorm(0.975) * pt$se
  names(theta) <- info$free_labels

  # Fit statistics -----------------------------------------------------------
  # Mplus counts the known-class proportions as free parameters (unless the
  # class-membership logits were fixed by the user) and includes their
  # contribution in the log-likelihood of KNOWNCLASS models.
  ll <- best$ll
  npar <- info$d
  K <- model$k
  int_free <- if (K > 1L) any(pt$free[pt$type == "intercept"]) else TRUE
  known_ll <- if (D$G > 1L) sum(D$n_g * log(D$n_g / D$N)) else 0
  known_par <- if (D$G > 1L && int_free) D$G - 1L else 0L
  # Mplus brings predictors of profile membership into the likelihood as
  # multivariate normal variables (saturated means and covariance matrix).
  pred_ll <- 0; pred_par <- 0L
  if (D$q > 0L) {
    Xc <- scale(D$X, scale = FALSE)
    S <- crossprod(Xc) / D$N
    pred_ll <- -D$N / 2 * (D$q * LOG2PI + as.numeric(determinant(S, logarithm = TRUE)$modulus) + D$q)
    pred_par <- as.integer(D$q + D$q * (D$q + 1) / 2)
  }
  ll_m <- ll + known_ll + pred_ll
  npar_m <- npar + known_par + pred_par
  fit <- if (mplus_convention) fit_indices(ll_m, npar_m, D$N) else fit_indices(ll, npar, D$N)

  tau <- es$tau
  colnames(tau) <- paste0("c", seq_len(K))
  cs <- class_summary(tau, es, D, K)
  reps <- sum(abs(res$final_ll - ll) < 1e-3)

  # Mplus-style estimation warnings -------------------------------------------
  warn <- estimation_warnings(converged = best$converged, iterations = best$iter, tol = tol, K = K,
                              n_final = length(res$final_ll), reps = reps, n_failed = res$n_failed,
                              singular = singular)
  if (length(warn))
    warning(paste(c(paste0("mglpa (", model_description(model), "):"), warn), collapse = "\n  "),
            call. = FALSE)

  structure(list(
    model = model,
    data = list(N = D$N, n_g = D$n_g, group_levels = D$group_levels,
                n_dropped = D$n_dropped, variables = D$vars,
                predictors = model$predictors),
    sample_stats = sample_statistics(D),
    D = D,
    info = info,
    data_frame = data_frame,
    fixed = fixed,
    coefficients = theta,
    params = pt,
    loglik = if (mplus_convention) ll_m else ll,
    loglik_mixture = ll,
    loglik_known_groups = known_ll,
    loglik_predictors = pred_ll,
    npar = if (mplus_convention) npar_m else npar,
    npar_free = npar,
    mplus_convention = mplus_convention,
    fit = fit,
    entropy = if (mplus_convention) cs$entropy_patterns else cs$entropy,
    entropy_profiles = cs$entropy,
    class_counts = cs,
    posterior = tau,
    class = cs$most_likely_class,
    starts = list(requested = starts, stiterations = stiterations, maxiter = maxiter, tol = tol,
                  final_ll = res$final_ll, converged = res$converged,
                  best_ll_replications = reps, stage1_ll = res$stage1_ll,
                  extra_iterations = best$extra_iterations, n_failed = res$n_failed),
    converged = best$converged,
    iterations = best$iter,
    warnings = warn,
    vcov = vc, vcov_ml = vc_ml, information_condition = cond, information_singular = singular,
    se_type = se,
    time = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ), class = "mglpa_fit")
}

# Means and standard deviations of the analysed variables, pooled and by group.
sample_statistics <- function(D) {
  pooled <- data.frame(variable = D$vars,
                       mean = colMeans(D$Y, na.rm = TRUE),
                       sd = apply(D$Y, 2, stats::sd, na.rm = TRUE),
                       n = colSums(D$O), stringsAsFactors = FALSE, row.names = NULL)
  by_group <- do.call(rbind, lapply(seq_len(D$G), function(g) {
    Yg <- D$Y[D$grp_rows[[g]], , drop = FALSE]
    data.frame(group = g, group_level = as.character(D$group_levels[g]), variable = D$vars,
               mean = colMeans(Yg, na.rm = TRUE), sd = apply(Yg, 2, stats::sd, na.rm = TRUE),
               n = colSums(!is.na(Yg)), stringsAsFactors = FALSE, row.names = NULL)
  }))
  list(pooled = pooled, by_group = by_group)
}

# ------------------------------------------------------------------------------
# Data preparation --------------------------------------------------------------
# ------------------------------------------------------------------------------

prepare_data <- function(data, model) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  vars <- c(model$indicators, model$outcomes)
  needed <- c(vars, model$predictors, model$group)
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars))
    stop("Variables not found in `data`: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  Y <- as.matrix(data[, vars, drop = FALSE])
  if (!is.numeric(Y)) stop("Indicators and outcomes must be numeric.", call. = FALSE)
  storage.mode(Y) <- "double"
  X <- NULL
  if (length(model$predictors)) {
    X <- as.matrix(data[, model$predictors, drop = FALSE])
    if (!is.numeric(X)) stop("Predictors must be numeric.", call. = FALSE)
    storage.mode(X) <- "double"
  }
  if (!is.null(model$group)) {
    gv <- data[[model$group]]
    lev <- model$group_levels
    if (is.null(lev)) lev <- sort(unique(gv[!is.na(gv)]))
    g <- match(gv, lev)
    if (any(is.na(g) & !is.na(gv)))
      stop("Some group values are not in `group_levels`.", call. = FALSE)
  } else {
    lev <- "all"
    g <- rep(1L, nrow(data))
  }
  n_ind <- length(model$indicators)
  keep <- rowSums(!is.na(Y[, seq_len(n_ind), drop = FALSE])) > 0 & !is.na(g)
  if (!is.null(X)) keep <- keep & rowSums(is.na(X)) == 0
  n_dropped <- sum(!keep)
  Y <- Y[keep, , drop = FALSE]
  g <- g[keep]
  if (!is.null(X)) X <- X[keep, , drop = FALSE]
  if (nrow(Y) < 2L) stop("Fewer than two usable observations.", call. = FALSE)
  G <- length(lev)
  grp_rows <- lapply(seq_len(G), function(i) which(g == i))
  n_g <- lengths(grp_rows)
  if (any(n_g == 0L))
    stop("Group(s) without observations: ", paste(lev[n_g == 0L], collapse = ", "), call. = FALSE)
  O <- !is.na(Y)
  Y0 <- Y
  Y0[!O] <- 0
  O1 <- O * 1
  Y0sq <- Y0 * Y0
  # Per-group slices used in every EM iteration
  grp <- lapply(seq_len(G), function(i) {
    r <- grp_rows[[i]]
    list(Y0 = Y0[r, , drop = FALSE], Y0sq = Y0sq[r, , drop = FALSE], O1 = O1[r, , drop = FALSE],
         X = if (is.null(X)) NULL else X[r, , drop = FALSE], n = length(r))
  })
  list(Y = Y, Y0 = Y0, Y0sq = Y0sq, O = O, O1 = O1, X = X, g = g, G = G,
       group_levels = lev, grp_rows = grp_rows, grp = grp, n_g = n_g, N = nrow(Y),
       p = ncol(Y), q = if (is.null(X)) 0L else ncol(X), n_dropped = n_dropped,
       vars = vars)
}

# ------------------------------------------------------------------------------
# Indexing information used by the E- and M-steps -------------------------------
# ------------------------------------------------------------------------------

make_info <- function(pt, D, model, engine = "auto") {
  K <- model$k; G <- D$G; p <- D$p; q <- D$q
  free_labels <- unique(pt$label[pt$free])
  pt$index <- NA_integer_
  pt$index[pt$free] <- match(pt$label[pt$free], free_labels)
  d <- length(free_labels)

  mean_rows <- pt$type == "mean"
  var_rows <- pt$type == "variance"
  mean_members <- split(pt$pos[mean_rows & pt$free], pt$label[mean_rows & pt$free])
  var_members <- split(pt$pos[var_rows & pt$free], pt$label[var_rows & pt$free])
  # Positions of the free parameters in the p x K x G arrays and integer ids
  # (1..L) of the label they belong to, for vectorised pooling with rowsum().
  agg <- function(members) {
    if (!length(members)) return(list(pos = integer(0), id = integer(0), n = 0L))
    pos <- unlist(members, use.names = FALSE)
    id <- rep(seq_along(members), lengths(members))
    o <- order(pos)
    list(pos = pos[o], id = id[o], n = length(members))
  }
  mean_agg <- agg(mean_members)
  var_agg <- agg(var_members)

  # First row of every free label: used to pack arrays back into theta
  first <- match(free_labels, pt$label)
  pack_type <- pt$type[first]
  pack_pos <- pt$pos[first]

  sv <- apply(D$Y, 2, stats::var, na.rm = TRUE)
  sv[!is.finite(sv) | sv <= 0] <- 1
  var_floor <- array(rep(sv, times = K * G) * 1e-4, c(p, K, G))

  logit_closed <- K == 1L   # no membership logits in a one-profile model
  Xd <- NULL; mask <- NULL; int_row <- NULL; slp_row <- NULL
  if (K > 1L) {
    int_rows <- pt$type == "intercept"
    logit_closed <- q == 0L && all(pt$free[int_rows])
    if (!logit_closed) {
      Xint <- if (model$equal_proportions) matrix(1, D$N, 1L) else
        vapply(seq_len(G), function(gg) as.numeric(D$g == gg), numeric(D$N))
      Xint <- matrix(Xint, nrow = D$N)
      Xslp <- NULL
      if (q > 0L) {
        Xslp <- if (model$predictor_effects == "equal") D$X else
          do.call(cbind, lapply(seq_len(G), function(gg) D$X * (D$g == gg)))
      }
      Xd <- cbind(Xint, Xslp)
      n_int <- ncol(Xint)
      int_row <- function(g) if (model$equal_proportions) 1L else as.integer(g)
      slp_row <- function(x, g) n_int + if (model$predictor_effects == "equal") as.integer(x)
        else (as.integer(g) - 1L) * q + as.integer(x)
      dd <- ncol(Xd)
      mask <- matrix(FALSE, dd, K - 1L)
      for (r in which(int_rows)) mask[int_row(pt$group[r]), pt$class[r]] <- pt$free[r]
      for (r in which(pt$type == "slope"))
        mask[slp_row(pt$var_index[r], pt$group[r]), pt$class[r]] <- pt$free[r]
    }
  }
  # Information for vectorised random starts: variable, sampling pool (0 = all
  # observations, g = group g) and type of every free parameter, plus pooled
  # and group-wise sample means and variances.
  n_groups_of <- vapply(free_labels, function(l) length(unique(pt$group[pt$label == l])), integer(1))
  pool <- ifelse(n_groups_of > 1L, 0L, pt$group[first])
  pool_mean <- cbind(colMeans(D$Y, na.rm = TRUE),
                     vapply(seq_len(G), function(g) colMeans(D$Y[D$grp_rows[[g]], , drop = FALSE], na.rm = TRUE), numeric(p)))
  pool_var <- cbind(sv, vapply(seq_len(G), function(g)
    apply(D$Y[D$grp_rows[[g]], , drop = FALSE], 2, stats::var, na.rm = TRUE), numeric(p)))
  pool_mean[!is.finite(pool_mean)] <- 0
  pool_var[!is.finite(pool_var) | pool_var <= 0] <- 1
  sel_m <- which(pack_type == "mean"); sel_v <- which(pack_type == "variance")
  sel_i <- which(pack_type == "intercept")
  start <- list(mean_idx = sel_m, mean_j = pt$var_index[first][sel_m], mean_pool = pool[sel_m],
                var_idx = sel_v, var_j = pt$var_index[first][sel_v], var_pool = pool[sel_v],
                int_idx = sel_i, int_k = pt$class[first][sel_i], int_g = pt$group[first][sel_i],
                pool_mean = pool_mean, pool_var = pool_var)

  # Compiled EM path: label ids (0 = fixed) over the p x K x G array positions
  use_cpp <- logit_closed && !identical(engine, "r") && !isTRUE(getOption("mglpa.engine") == "r")
  cpp <- NULL
  if (use_cpp) {
    ids <- function(members) {
      id <- integer(p * K * G)
      for (i in seq_along(members)) id[members[[i]]] <- i
      id
    }
    cpp <- list(mean_id = ids(mean_members), var_id = ids(var_members),
                n_mean = length(mean_members), n_var = length(var_members))
  }
  list(pt = pt, free_labels = free_labels, d = d, mean_members = mean_members,
       var_members = var_members, mean_agg = mean_agg, var_agg = var_agg,
       pack_type = pack_type, pack_pos = pack_pos,
       var_floor = var_floor, logit_closed = logit_closed, Xd = Xd,
       logit_mask = mask, int_row = int_row, slp_row = slp_row,
       use_cpp = use_cpp, cpp = cpp, start = start,
       K = K, G = G, p = p, q = q)
}

unpack_theta <- function(theta, info) {
  pt <- info$pt
  vals <- pt$value
  vals[pt$free] <- theta[pt$index[pt$free]]
  vals_to_par(vals, info)
}

vals_to_par <- function(vals, info) {
  pt <- info$pt
  K <- info$K; G <- info$G; p <- info$p; q <- info$q
  list(
    mu = array(vals[pt$type == "mean"], c(p, K, G)),
    sig2 = array(vals[pt$type == "variance"], c(p, K, G)),
    alpha = if (K > 1L) array(vals[pt$type == "intercept"], c(K - 1L, G)) else NULL,
    beta = if (K > 1L && q > 0L) array(vals[pt$type == "slope"], c(q, K - 1L, G)) else NULL
  )
}

pack_theta <- function(par, info) {
  theta <- numeric(info$d)
  for (type in c("mean", "variance", "intercept", "slope")) {
    sel <- info$pack_type == type
    if (!any(sel)) next
    arr <- switch(type, mean = par$mu, variance = par$sig2, intercept = par$alpha,
                  slope = par$beta)
    theta[sel] <- arr[info$pack_pos[sel]]
  }
  theta
}

# ------------------------------------------------------------------------------
# E-step ----------------------------------------------------------------------
# ------------------------------------------------------------------------------

estep <- function(par, D, info) {
  N <- D$N; K <- info$K; G <- info$G; p <- info$p
  logf <- matrix(0, N, K)
  logpi <- matrix(0, N, K)
  for (gg in seq_len(G)) {
    rows <- D$grp_rows[[gg]]
    S <- D$grp[[gg]]
    ng <- S$n
    # Sum over observed indicators of log N(y_j; mu_jk, s2_jk), for all
    # profiles at once: three matrix products of (n x p) by (p x K).
    mu <- matrix(par$mu[, , gg], nrow = p)
    s2 <- matrix(par$sig2[, , gg], nrow = p)
    inv <- 1 / s2
    const <- LOG2PI + log(s2) + mu * mu * inv
    logf[rows, ] <- -0.5 * (S$O1 %*% const - 2 * (S$Y0 %*% (mu * inv)) + S$Y0sq %*% inv)
    if (K > 1L) {
      eta <- matrix(0, ng, K)
      eta[, -K] <- rep(par$alpha[, gg], each = ng)
      if (!is.null(par$beta)) {
        eta[, -K] <- eta[, -K] + S$X %*% matrix(par$beta[, , gg], nrow = info$q, ncol = K - 1L)
      }
      mx <- eta[cbind(seq_len(ng), max.col(eta, ties.method = "first"))]
      lse <- mx + log(rowSums(exp(eta - mx)))
      logpi[rows, ] <- eta - lse
    }
  }
  lp <- logf + logpi
  if (K > 1L) {
    mx <- lp[cbind(seq_len(N), max.col(lp, ties.method = "first"))]
    ll_i <- mx + log(rowSums(exp(lp - mx)))
    tau <- exp(lp - ll_i)
  } else {
    ll_i <- lp[, 1L]
    tau <- matrix(1, N, 1L)
  }
  list(tau = tau, ll_i = ll_i, ll = sum(ll_i), logpi = logpi, logf = logf)
}

# ------------------------------------------------------------------------------
# M-step (conditional maximisation with pooled sufficient statistics) ---------
# ------------------------------------------------------------------------------

mstep <- function(tau, par, D, info, model) {
  p <- info$p; K <- info$K; G <- info$G
  A <- B <- C <- array(0, c(p, K, G))
  for (gg in seq_len(G)) {
    rows <- D$grp_rows[[gg]]
    S <- D$grp[[gg]]
    tg <- tau[rows, , drop = FALSE]
    A[, , gg] <- crossprod(S$Y0, tg)
    B[, , gg] <- crossprod(S$O1, tg)
    C[, , gg] <- crossprod(S$Y0sq, tg)
  }
  mu <- par$mu
  s2 <- par$sig2
  # Means: precision-weighted pooling over the positions sharing a label
  ma <- info$mean_agg
  if (ma$n > 0L) {
    w <- 1 / s2[ma$pos]
    num <- rowsum(A[ma$pos] * w, ma$id, reorder = TRUE)[, 1L]
    den <- rowsum(B[ma$pos] * w, ma$id, reorder = TRUE)[, 1L]
    ok <- den > 0
    val <- mu[ma$pos]
    val[ok[ma$id]] <- (num / den)[ma$id][ok[ma$id]]
    mu[ma$pos] <- val
  }
  # Variances: pooled sums of squares around the updated means
  va <- info$var_agg
  if (va$n > 0L) {
    SS <- C - 2 * mu * A + mu * mu * B
    num <- rowsum(SS[va$pos], va$id, reorder = TRUE)[, 1L]
    den <- rowsum(B[va$pos], va$id, reorder = TRUE)[, 1L]
    ok <- den > 0
    val <- s2[va$pos]
    new <- pmax(num / den, info$var_floor[va$pos][match(seq_len(va$n), va$id)])
    val[ok[va$id]] <- new[va$id][ok[va$id]]
    s2[va$pos] <- val
  }
  alpha <- par$alpha
  beta <- par$beta
  if (K > 1L) {
    if (info$logit_closed) {
      if (model$equal_proportions) {
        pk <- pmax(colSums(tau) / D$N, 1e-12)
        alpha[] <- log(pk[-K] / pk[K])
      } else {
        for (gg in seq_len(G)) {
          rows <- D$grp_rows[[gg]]
          pk <- pmax(colSums(tau[rows, , drop = FALSE]) / length(rows), 1e-12)
          alpha[, gg] <- log(pk[-K] / pk[K])
        }
      }
    } else {
      th <- logit_theta_from_par(alpha, beta, info)
      th <- logit_newton(tau, info$Xd, th, info$logit_mask)
      pb <- logit_par_from_theta(th, info)
      alpha <- pb$alpha
      beta <- pb$beta
    }
  }
  list(mu = mu, sig2 = s2, alpha = alpha, beta = beta)
}

# ------------------------------------------------------------------------------
# EM loop ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

em_iterate <- function(par, D, info, model, maxiter, tol, check = TRUE) {
  if (isTRUE(info$use_cpp)) return(em_iterate_cpp(par, D, info, model, maxiter, tol, check))
  ll_old <- -Inf
  converged <- FALSE
  es <- NULL
  iter <- 0L
  for (it in seq_len(maxiter)) {
    es <- estep(par, D, info)
    iter <- it
    if (!is.finite(es$ll)) break
    if (check && abs(es$ll - ll_old) < tol) {
      converged <- TRUE
      break
    }
    ll_old <- es$ll
    par <- mstep(es$tau, par, D, info, model)
  }
  if (is.null(es)) es <- estep(par, D, info)
  # `es` was computed at `par` before the last M-step; recompute so that the
  # returned parameters and posteriors are consistent.
  if (!converged || it == maxiter) es <- estep(par, D, info)
  list(par = par, es = es, ll = es$ll, iter = iter, converged = converged)
}

# A PSOCK cluster whose workers have the package namespace loaded, or NULL
# (with a warning) when the package is not installed where the workers can
# find it, e.g. during development with pkgload::load_all().
start_cluster <- function(cores) {
  cl <- tryCatch(parallel::makeCluster(cores), error = function(e) NULL)
  if (is.null(cl)) {
    warning("Could not start a parallel cluster; estimating sequentially.", call. = FALSE)
    return(NULL)
  }
  lp <- .libPaths()
  ok <- tryCatch({
    parallel::clusterCall(cl, function(lp) { .libPaths(lp); invisible(NULL) }, lp)
    unlist(parallel::clusterEvalQ(cl, requireNamespace("mglpa", quietly = TRUE)))
  }, error = function(e) FALSE)
  if (!all(ok)) {
    parallel::stopCluster(cl)
    warning("Parallel estimation needs the installed mglpa package on the worker processes; ",
            "estimating sequentially.", call. = FALSE)
    return(NULL)
  }
  cl
}

# Same iterations in compiled code (see src/em.cpp). Profile proportions are
# passed as probabilities and converted back to logits relative to the last
# profile. The posterior probabilities (`es`) are computed on demand by the
# caller, so `es` is NULL here.
em_iterate_cpp <- function(par, D, info, model, maxiter, tol, check = TRUE) {
  K <- info$K; G <- info$G
  pi <- matrix(1, K, G)
  if (K > 1L) for (g in seq_len(G)) {
    a <- c(par$alpha[, g], 0)
    e <- exp(a - max(a))
    pi[, g] <- e / sum(e)
  }
  res <- em_core_cpp(D$Y0, D$O1, D$g, G, as.numeric(par$mu), as.numeric(par$sig2), pi,
                     isTRUE(model$equal_proportions), info$cpp$mean_id, info$cpp$var_id,
                     info$cpp$n_mean, info$cpp$n_var, as.numeric(info$var_floor),
                     as.integer(max(maxiter, 1L)), tol, check)
  alpha <- NULL
  if (K > 1L) {
    pk <- res$pi
    alpha <- log(pk[-K, , drop = FALSE] / rep(pk[K, ], each = K - 1L))
  }
  par_new <- list(mu = array(res$mu, dim(par$mu)), sig2 = array(res$s2, dim(par$sig2)),
                  alpha = alpha, beta = NULL)
  list(par = par_new, es = NULL, ll = res$ll, iter = res$iter, converged = res$converged)
}

# ------------------------------------------------------------------------------
# Random starts ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# Random starting values (Mplus-like perturbation): every free profile mean is
# set to the value of a randomly drawn observation of its group(s), variances
# to the sample variance times a random factor in [0.5, 1.5], proportions to
# random draws, and regression slopes to zero. Returns a value for every row
# of the parameter table (fixed rows keep their fixed value).
random_start <- function(pt, D, info, model) {
  si <- info$start
  K <- info$K; G <- info$G
  lv <- numeric(info$d)
  if (length(si$mean_idx)) {
    rows <- integer(length(si$mean_idx))
    for (pl in unique(si$mean_pool)) {
      sel <- si$mean_pool == pl
      cand <- if (pl == 0L) seq_len(D$N) else D$grp_rows[[pl]]
      rows[sel] <- cand[sample.int(length(cand), sum(sel), replace = TRUE)]
    }
    yv <- D$Y[cbind(rows, si$mean_j)]
    na <- is.na(yv)
    if (any(na)) yv[na] <- si$pool_mean[cbind(si$mean_j[na], si$mean_pool[na] + 1L)]
    lv[si$mean_idx] <- yv
  }
  if (length(si$var_idx)) {
    lv[si$var_idx] <- si$pool_var[cbind(si$var_j, si$var_pool + 1L)] *
      stats::runif(length(si$var_idx), 0.5, 1.5)
  }
  if (length(si$int_idx)) {
    alpha0 <- matrix(0, K - 1L, G)
    if (model$equal_proportions) {
      w <- stats::runif(K, 0.2, 1)
      pk <- w / sum(w)
      alpha0[] <- log(pk[-K] / pk[K])
    } else {
      for (gg in seq_len(G)) {
        w <- stats::runif(K, 0.2, 1)
        pk <- w / sum(w)
        alpha0[, gg] <- log(pk[-K] / pk[K])
      }
    }
    lv[si$int_idx] <- alpha0[cbind(si$int_k, si$int_g)]
  }
  vals <- pt$value
  vals[pt$free] <- lv[pt$index[pt$free]]
  vals
}

run_random_starts <- function(pt, D, info, model, starts, stiterations, maxiter,
                              tol, user_vals, verbose, cores = 1L) {
  n1 <- as.integer(starts[1L]); n2 <- as.integer(starts[2L])
  inits <- list()
  if (!is.null(user_vals)) inits[[1L]] <- user_vals
  for (s in seq_len(n1)) inits[[length(inits) + 1L]] <- random_start(pt, D, info, model)
  if (!length(inits)) inits[[1L]] <- random_start(pt, D, info, model)
  n_init <- length(inits)
  cl <- NULL
  cores <- max(1L, as.integer(cores))
  if (cores > 1L && n_init > 1L) {
    cl <- start_cluster(cores)
    if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  run <- function(X, FUN) if (is.null(cl)) lapply(X, FUN) else parallel::parLapply(cl, X, FUN)
  if (verbose) message("Initial stage: ", n_init, " starts x ", stiterations, " iterations")
  stage1_fun <- function(vals) {
    r <- em_iterate(vals_to_par(vals, info), D, info, model, maxiter = max(stiterations, 1L),
                    tol = tol, check = FALSE)
    list(par = r$par, ll = r$ll)
  }
  s1 <- run(inits, stage1_fun)
  stage1 <- lapply(s1, `[[`, "par")
  ll1 <- vapply(s1, `[[`, numeric(1), "ll")
  ll1[!is.finite(ll1)] <- -Inf
  ord <- order(ll1, decreasing = TRUE)
  keep <- ord[seq_len(min(max(n2, 1L), n_init))]
  if (verbose) message("Final stage: ", length(keep), " optimisations")
  stage2 <- run(stage1[keep], function(par) em_iterate(par, D, info, model, maxiter, tol))
  ll2 <- vapply(stage2, `[[`, numeric(1), "ll")
  ll2[!is.finite(ll2)] <- -Inf
  if (all(ll2 == -Inf))
    stop("The model could not be estimated: every random start produced a non-finite ",
         "log-likelihood. Check the data (constant or extreme indicators) and the model.",
         call. = FALSE)
  best_i <- which.max(ll2)
  best <- stage2[[best_i]]
  best$extra_iterations <- 0L
  if (!best$converged) {
    # Give the best solution more iterations before declaring non-convergence
    more <- em_iterate(best$par, D, info, model, maxiter * 4L, tol)
    best$extra_iterations <- more$iter
    best$iter <- best$iter + more$iter
    best$par <- more$par; best$es <- more$es; best$ll <- more$ll
    best$converged <- more$converged
    ll2[best_i] <- best$ll
  }
  if (is.null(best$es)) {
    best$es <- estep(best$par, D, info)
    best$ll <- best$es$ll
  }
  conv <- vapply(stage2, `[[`, logical(1), "converged")
  conv[best_i] <- best$converged
  ord2 <- order(ll2, decreasing = TRUE)
  list(best = best, final_ll = ll2[ord2], converged = conv[ord2], stage1_ll = ll1,
       n_failed = sum(ll2 == -Inf))
}
