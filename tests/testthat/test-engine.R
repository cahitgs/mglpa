test_that("a one-profile model gives the closed-form normal estimates", {
  fit <- fit_lpa(wm1, mglpa_model(1, ind), starts = c(1, 1), seed = 1)
  Y <- as.matrix(wm1[, ind])
  n <- nrow(Y)
  expect_equal(unname(coef(fit)[paste0("m_", ind, "_c1_g1")]), unname(colMeans(Y)), tolerance = 1e-6)
  expect_equal(unname(coef(fit)[paste0("v_", ind, "_g1")]),
               unname(apply(Y, 2, function(x) mean((x - mean(x))^2))), tolerance = 1e-5)
  ll <- sum(vapply(seq_along(ind), function(j) {
    x <- Y[, j]; v <- mean((x - mean(x))^2)
    sum(dnorm(x, mean(x), sqrt(v), log = TRUE))
  }, numeric(1)))
  expect_equal(fit$loglik, ll, tolerance = 1e-6)
  expect_equal(fit$npar, 10)
  expect_equal(fit$loglik, ref_value("enumeration", "ll", group = 1, k = 1), tolerance = 1e-5)
  expect_true(is.na(fit$entropy))
})

test_that("the three-profile solution of group 1 reproduces Mplus", {
  fit <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(30, 5), seed = 11)
  expect_equal(fit$loglik, ref_value("enumeration", "ll", group = 1, k = 3), tolerance = 1e-5)
  expect_equal(fit$npar, 22)
  expect_equal(fit$fit[["BIC"]], ref_value("enumeration", "bic", group = 1, k = 3), tolerance = 1e-5)
  expect_equal(fit$fit[["CAIC"]], ref_value("enumeration", "caic", group = 1, k = 3), tolerance = 1e-5)
  expect_equal(fit$fit[["ABIC"]], ref_value("enumeration", "abic", group = 1, k = 3), tolerance = 1e-5)
  expect_equal(fit$entropy, 0.974, tolerance = 2e-3)
  expect_equal(sort(fit$class_counts$most_likely$proportion), c(0.208, 0.356, 0.436), tolerance = 1e-6)
  expect_true(fit$converged)
  expect_equal(fit$starts$best_ll_replications, 5)

  # MLR standard errors (Mplus: ESTIMATOR = MLR) for a few parameters
  pm <- parameters(fit, "mean")
  lo <- pm[pm$variable == "verbal_wm" & abs(pm$est + 1.152) < 0.01, ]
  expect_equal(lo$se, 0.051, tolerance = 0.05)
  hi <- pm[pm$variable == "verbal_wm" & abs(pm$est - 1.147) < 0.01, ]
  expect_equal(hi$se, 0.030, tolerance = 0.05)
  pv <- parameters(fit, "variance")
  expect_equal(pv$se[pv$variable == "verbal_wm"][1], 0.014, tolerance = 0.08)
  expect_equal(pv$est[pv$variable == "verbal_wm"][1], 0.208, tolerance = 5e-3)

  # The analytic gradient vanishes at the maximum
  g <- mglpa:::loglik_gradient(coef(fit), fit$D, fit$info)
  expect_lt(max(abs(g)), 1e-2)

  # ML standard errors are also available
  fit_ml <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(10, 2), seed = 11, se = "ml")
  expect_true(all(is.finite(parameters(fit_ml)$se)))
})

test_that("missing indicator values are handled by FIML", {
  set.seed(5)
  w <- wm1
  w$verbal_wm[sample(nrow(w), 50)] <- NA
  w$updating[sample(nrow(w), 30)] <- NA
  fit <- fit_lpa(w, mglpa_model(2, ind), starts = c(10, 2), seed = 3)
  expect_true(is.finite(fit$loglik))
  expect_true(fit$converged)
  expect_equal(fit$data$N, nrow(w))
  expect_equal(nrow(fit$posterior), nrow(w))
  expect_true(all(abs(rowSums(fit$posterior) - 1) < 1e-10))
})

test_that("class-varying variances are estimated per profile", {
  fit <- fit_lpa(wm1, mglpa_model(2, ind, variances = "class_varying"), starts = c(10, 2), seed = 2)
  expect_equal(fit$npar, 5 * 2 + 5 * 2 + 1)
  v <- parameters(fit, "variance")
  expect_false(isTRUE(all.equal(v$est[v$class == 1], v$est[v$class == 2])))
})

test_that("relabelling profiles keeps the likelihood and reorders the estimates", {
  fit <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(20, 4), seed = 4)
  f_level <- reorder_profiles(fit, "level")
  pm <- profile_means(f_level)
  lev <- vapply(1:3, function(k) mean(pm$z[pm$class == k]), numeric(1))
  expect_true(all(diff(lev) > 0))
  expect_equal(f_level$loglik, fit$loglik, tolerance = 1e-6)
  f_size <- reorder_profiles(fit, "size")
  expect_true(all(diff(f_size$class_counts$by_group_model[1, ]) <= 0))
  f_perm <- reorder_profiles(f_level, c(3, 1, 2))
  expect_equal(unname(coef(f_perm)["m_verbal_wm_c1_g1"]), unname(coef(f_level)["m_verbal_wm_c3_g1"]))
  expect_equal(f_perm$loglik, fit$loglik, tolerance = 1e-6)
  expect_equal(f_perm$class_counts$by_group_model[1, 1], f_level$class_counts$by_group_model[1, 3],
               tolerance = 1e-6)
})

test_that("fixed parameters are held constant and not counted", {
  fit <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(20, 4), seed = 4)
  fx <- fixed_from_fit(fit, c("mean", "variance"))
  fit2 <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(2, 1), fixed = fx, start_values = fit)
  expect_equal(fit2$npar, 2)
  expect_equal(fit2$loglik, fit$loglik, tolerance = 1e-6)
  expect_equal(unname(coef(fit2)), unname(coef(fit)[names(coef(fit2))]), tolerance = 1e-4)
  expect_false(any(parameters(fit2)$free[parameters(fit2)$type == "mean"]))
})

test_that("posterior probabilities and summaries are consistent", {
  fit <- fit_lpa(wm1, mglpa_model(2, ind), starts = c(10, 2), seed = 6)
  pp <- posterior_probabilities(fit)
  expect_equal(ncol(pp), 3)
  expect_equal(pp$most_likely, max.col(as.matrix(pp[, 1:2])))
  expect_equal(sum(fit$class_counts$estimated$count), nrow(wm1))
  expect_equal(sum(fit$class_counts$most_likely$count), nrow(wm1))
  expect_s3_class(logLik(fit), "logLik")
  expect_equal(nobs(fit), nrow(wm1))
  expect_output(print(fit), "Log-likelihood")
  expect_output(print(summary(fit)), "Parameters")
})
