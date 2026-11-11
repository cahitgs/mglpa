#' Vuong-Lo-Mendell-Rubin likelihood ratio tests for the number of profiles
#'
#' Compares a k-profile solution with the (k - 1)-profile solution estimated
#' on the same data with the Vuong-Lo-Mendell-Rubin likelihood ratio test
#' (VLMR; Lo, Mendell and Rubin, 2001) and its adjusted version (aLMR), the
#' two tests reported in the Mplus `TECH11` output.
#'
#' @details
#' The test statistic is \eqn{2(LL_k - LL_{k-1})}. Following Vuong (1989) and
#' Lo et al. (2001), its null distribution is approximated by a weighted sum
#' of independent chi-square variables with one degree of freedom whose
#' weights are the eigenvalues of the matrix
#' \deqn{W = \left(\begin{array}{cc} C_k V_k & C_{k,k-1} V_{k-1} \\
#'   -C_{k-1,k} V_k & -C_{k-1} V_{k-1} \end{array}\right),}
#' where the \eqn{C} matrices are cross-products of the per-observation score
#' vectors of the two models at their maximum likelihood estimates and the
#' \eqn{V} matrices are covariance matrices of the parameter estimates.
#' Vuong's derivation uses the inverse observed information matrices
#' (`variant = "vuong"`, the implementation of LatentGOLD). Mplus uses the
#' robust sandwich covariance matrices instead (Vermunt, 2024);
#' `variant = "mplus"`, the default, reproduces the `TECH11` output of Mplus
#' (statistic, mean, standard deviation and p-values of the VLMR and aLMR
#' tests). The mean of the reference distribution is the sum of the
#' eigenvalues, its standard deviation the square root of twice the sum of
#' their squares, and p-values are computed with the method of Imhof (1961).
#' The adjusted statistic of Lo et al. (2001) divides the statistic by
#' \eqn{1 + 1 / (d \log N)}, where \eqn{d} is the difference in the number of
#' free parameters and \eqn{N} the sample size, and uses the same reference
#' distribution.
#'
#' As in Mplus, the tests are available for models estimated in a single
#' group without predictors; the two fits must differ only in the number of
#' profiles. Jeffries (2003) showed that the derivation of Lo et al. is not
#' rigorous, so these p-values are best used together with the bootstrap
#' likelihood ratio test ([blrt()]) and the information criteria.
#'
#' @param fit_k,fit_km1 `mglpa_fit` objects with k and k - 1 profiles
#'   estimated on the same single-group data with the same variance
#'   specification.
#' @param variant `"mplus"` (sandwich covariance matrices, as in the Mplus
#'   `TECH11` output) or `"vuong"` (inverse information matrices, Vuong's
#'   original matrix).
#' @return An object of class `mglpa_vlmr`: a list with the number of
#'   profiles `k`, the `statistic`, the difference in the number of parameters
#'   `df`, the sample size `N`, the two log-likelihoods, the `mean`, `sd` and
#'   `p` of the VLMR test, the `adjusted` statistic and `adjusted_p` of the
#'   aLMR test, and the `eigenvalues` of the weight matrix.
#' @references
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in
#' normal variables. *Biometrika, 48*, 419-426.
#'
#' Jeffries, N. O. (2003). A note on 'Testing the number of components in a
#' normal mixture'. *Biometrika, 90*, 991-994.
#'
#' Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
#' components in a normal mixture. *Biometrika, 88*, 767-778.
#'
#' Vermunt, J. K. (2024). The Vuong-Lo-Mendell-Rubin test for latent class and
#' latent profile analysis: A note on the different implementations in Mplus
#' and LatentGOLD. *Methodology, 20*(1), Article e12467.
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. *Econometrica, 57*, 307-333.
#' @examples
#' \donttest{
#' data(mglpa_wm)
#' ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
#' d1 <- mglpa_wm[mglpa_wm$group == 1, ]
#' f2 <- fit_lpa(d1, mglpa_model(2, ind), starts = c(20, 4), seed = 1)
#' f3 <- fit_lpa(d1, mglpa_model(3, ind), starts = c(20, 4), seed = 1)
#' vlmr_test(f3, f2)
#' }
#' @export
vlmr_test <- function(fit_k, fit_km1, variant = c("mplus", "vuong")) {
  variant <- match.arg(variant)
  stopifnot(inherits(fit_k, "mglpa_fit"), inherits(fit_km1, "mglpa_fit"))
  if (fit_k$model$k != fit_km1$model$k + 1L)
    stop("`fit_k` must have one more profile than `fit_km1`.")
  if (fit_k$D$G > 1L || fit_km1$D$G > 1L || length(fit_k$model$predictors) ||
      length(fit_km1$model$predictors))
    stop("vlmr_test() is available for single-group models without predictors ",
         "(as TECH11 in Mplus).")
  if (fit_k$data$N != fit_km1$data$N || !identical(fit_k$D$vars, fit_km1$D$vars))
    stop("The two fits must be estimated on the same data.")
  N <- fit_k$data$N
  a <- fit_scores_information(fit_k)
  b <- fit_scores_information(fit_km1)
  Ca <- crossprod(a$scores)
  Cb <- crossprod(b$scores)
  Cab <- crossprod(a$scores, b$scores)
  Va <- safe_inverse(a$information)
  Vb <- safe_inverse(b$information)
  singular <- isTRUE(attr(Va, "singular")) || isTRUE(attr(Vb, "singular"))
  attr(Va, "singular") <- NULL
  attr(Vb, "singular") <- NULL
  if (variant == "mplus") {
    Va <- Va %*% Ca %*% Va
    Vb <- Vb %*% Cb %*% Vb
  }
  W <- rbind(cbind(Ca %*% Va, Cab %*% Vb),
             cbind(-t(Cab) %*% Va, -Cb %*% Vb))
  lambda <- Re(eigen(W, only.values = TRUE)$values)
  stat <- 2 * (fit_k$loglik_mixture - fit_km1$loglik_mixture)
  df <- length(fit_k$coefficients) - length(fit_km1$coefficients)
  adjusted <- stat / (1 + 1 / (df * log(N)))
  if (singular)
    warning("The information matrix of one of the models is singular or ill-conditioned; ",
            "the VLMR test may not be trustworthy.", call. = FALSE)
  structure(list(k = fit_k$model$k, statistic = stat, df = df, N = N,
                 loglik_k = fit_k$loglik_mixture, loglik_km1 = fit_km1$loglik_mixture,
                 mean = sum(lambda), sd = sqrt(2 * sum(lambda^2)),
                 p = imhof_p(stat, lambda),
                 adjusted = adjusted, adjusted_p = imhof_p(adjusted, lambda),
                 eigenvalues = lambda, variant = variant, singular = singular),
            class = "mglpa_vlmr")
}

#' @export
print.mglpa_vlmr <- function(x, digits = 3, ...) {
  cat(sprintf("Vuong-Lo-Mendell-Rubin likelihood ratio test for %d (H0) versus %d profiles\n",
              x$k - 1L, x$k))
  cat(sprintf("  H0 log-likelihood = %.*f\n", digits, x$loglik_km1))
  cat(sprintf("  2 x log-likelihood difference = %.*f, difference in parameters = %d\n",
              digits, x$statistic, x$df))
  cat(sprintf("  Mean = %.*f, SD = %.*f, p = %.4f\n", digits, x$mean, digits, x$sd, x$p))
  cat(sprintf("Lo-Mendell-Rubin adjusted LRT: value = %.*f, p = %.4f\n", digits, x$adjusted,
              x$adjusted_p))
  if (x$variant == "vuong")
    cat("  (Vuong's matrix with inverse information matrices; Mplus TECH11 uses sandwich",
        "covariance matrices)\n")
  invisible(x)
}

# Per-observation scores and observed information matrix of a fit at its
# maximum likelihood estimates.
fit_scores_information <- function(fit) {
  info <- fit$info
  inf <- observed_information(fit$coefficients, info$pt, fit$D, info, fit$model)
  list(scores = inf$scores, information = inf$information, singular = inf$singular)
}

# P(sum_j lambda_j Z_j^2 > x) for independent standard normal Z_j, by the
# numerical integration of Imhof (1961); falls back to simulation if the
# integration fails.
imhof_p <- function(x, lambda) {
  lambda <- lambda[is.finite(lambda) & abs(lambda) > 1e-10]
  if (!length(lambda) || !is.finite(x)) return(NA_real_)
  m <- sum(lambda)
  s <- sqrt(2 * sum(lambda^2))
  if (x > m + 60 * s) return(0)
  if (x < m - 60 * s) return(1)
  f <- function(u) {
    th <- 0.5 * colSums(atan(outer(lambda, u))) - 0.5 * x * u
    rho <- exp(0.25 * colSums(log1p(outer(lambda^2, u^2))))
    sin(th) / (u * rho)
  }
  val <- tryCatch(stats::integrate(f, 0, Inf, subdivisions = 2000L, rel.tol = 1e-9,
                                   abs.tol = 1e-12, stop.on.error = FALSE)$value,
                  error = function(e) NA_real_)
  if (is.finite(val)) return(min(max(0.5 + val / pi, 0), 1))
  # Simulation fallback (does not disturb the user's random number stream)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(1L)
  q <- numeric(0)
  for (i in seq_len(20L))
    q <- c(q, as.vector(matrix(stats::rnorm(5000L * length(lambda))^2, 5000L) %*% lambda))
  mean(q >= x)
}

#' Add Vuong-Lo-Mendell-Rubin tests to a class enumeration
#'
#' Runs [vlmr_test()] for every consecutive pair of solutions in each group
#' and adds the columns `LRT` (twice the log-likelihood difference), `VLMR_p`
#' and `aLMR_p` to the enumeration table. [enumerate_profiles()] calls this
#' function by default (`lrt = TRUE`).
#'
#' @param enum An `mglpa_enumeration` object.
#' @inheritParams vlmr_test
#' @return The updated `mglpa_enumeration` object, with the test objects in
#'   the element `vlmr`.
#' @export
enumeration_vlmr <- function(enum, variant = c("mplus", "vuong")) {
  stopifnot(inherits(enum, "mglpa_enumeration"))
  variant <- match.arg(variant)
  tab <- enum$table
  tab$LRT <- tab$VLMR_p <- tab$aLMR_p <- NA_real_
  enum$vlmr <- list()
  for (g in names(enum$fits)) {
    ks <- as.integer(names(enum$fits[[g]]))
    for (k in ks[-1]) {
      if (!(k - 1L) %in% ks) next
      res <- vlmr_test(enum$fits[[g]][[as.character(k)]], enum$fits[[g]][[as.character(k - 1L)]],
                       variant = variant)
      enum$vlmr[[paste(g, k)]] <- res
      i <- tab$group == g & tab$k == k
      tab$LRT[i] <- res$statistic
      tab$VLMR_p[i] <- res$p
      tab$aLMR_p[i] <- res$adjusted_p
    }
  }
  enum$table <- order_enum_columns(tab)
  enum
}

# Canonical column order of the enumeration table
order_enum_columns <- function(tab) {
  canonical <- c("model", "group", "k", "LL", "npar", "AIC", "BIC", "CAIC", "ABIC", "entropy",
                 "LRT", "VLMR_p", "aLMR_p", "BLRT_p", "min_proportion", "replications",
                 "n_final", "converged", "notes")
  tab[, c(intersect(canonical, names(tab)), setdiff(names(tab), canonical)), drop = FALSE]
}
