test_that("the Imhof integral reproduces chi-square probabilities", {
  expect_equal(mglpa:::imhof_p(7.81, rep(1, 3)), pchisq(7.81, 3, lower.tail = FALSE), tolerance = 1e-6)
  expect_equal(mglpa:::imhof_p(3, c(2, 2)), pchisq(1.5, 2, lower.tail = FALSE), tolerance = 1e-6)
  expect_equal(mglpa:::imhof_p(0.5, c(1, 1, 1, 1, 1)), pchisq(0.5, 5, lower.tail = FALSE), tolerance = 1e-6)
  expect_equal(mglpa:::imhof_p(1e6, rep(1, 3)), 0)
  expect_equal(mglpa:::imhof_p(-1e6, rep(1, 3)), 1)
  expect_true(is.na(mglpa:::imhof_p(1, numeric(0))))
})

test_that("VLMR and aLMR tests reproduce the Mplus TECH11 output", {
  skip_on_cran()
  t11 <- mplus_t11[mplus_t11$mplus_version == "8.3", ]
  check <- function(r, m) {
    expect_equal(r$loglik_km1, m$h0_ll, tolerance = 1e-5)
    expect_equal(r$statistic, m$lr, tolerance = 1e-4)
    expect_equal(r$df, m$df)
    expect_equal(r$mean, m$mean, tolerance = 2e-3)
    expect_equal(r$sd, m$sd, tolerance = 2e-3)
    expect_equal(r$p, m$vlmr_p, tolerance = 2e-3)
    expect_equal(r$adjusted, m$almr_value, tolerance = 1e-4)
    expect_equal(r$adjusted_p, m$almr_p, tolerance = 2e-3)
  }
  for (g in 1:2) {
    d <- wm[wm$group == g, ]
    fits <- lapply(1:4, function(k) fit_lpa(d, mglpa_model(k, ind), starts = c(300, 30),
                                             seed = 20 + k, se = "none", verbose = FALSE))
    for (k in 2:4) {
      r <- vlmr_test(fits[[k]], fits[[k - 1]])
      expect_s3_class(r, "mglpa_vlmr")
      check(r, t11[t11$group == g & t11$k == k, ])
      expect_output(print(r), "Vuong-Lo-Mendell-Rubin")
    }
  }
  # Vuong's original matrix (inverse information) gives a different reference distribution
  r_v <- vlmr_test(fits[[4]], fits[[3]], variant = "vuong")
  r_m <- vlmr_test(fits[[4]], fits[[3]])
  expect_equal(r_v$statistic, r_m$statistic)
  expect_false(isTRUE(all.equal(r_v$mean, r_m$mean)))
  expect_output(print(r_v), "inverse information")
  expect_error(vlmr_test(fits[[3]], fits[[4]]), "one more profile")
})

test_that("the VLMR test reproduces Mplus 9 with class-varying variances", {
  skip_on_cran()
  m <- mplus_t11[mplus_t11$mplus_version == "9.1.1" & mplus_t11$variances == "class_varying", ]
  skip_if(nrow(m) == 0)
  f2 <- fit_lpa(wm1, mglpa_model(2, ind, variances = "class_varying"), starts = c(300, 30),
                stiterations = 20, seed = 31, se = "none", verbose = FALSE)
  f3 <- fit_lpa(wm1, mglpa_model(3, ind, variances = "class_varying"), starts = c(300, 30),
                stiterations = 20, seed = 32, se = "none", verbose = FALSE)
  r <- vlmr_test(f3, f2)
  expect_equal(r$loglik_km1, m$h0_ll, tolerance = 1e-5)
  expect_equal(r$statistic, m$lr, tolerance = 1e-4)
  expect_equal(r$df, 11)
  expect_equal(r$mean, m$mean, tolerance = 2e-3)
  expect_equal(r$sd, m$sd, tolerance = 2e-3)
  expect_equal(r$adjusted, m$almr_value, tolerance = 1e-4)
  expect_equal(r$p, 0, tolerance = 1e-4)
})

test_that("the tests are restricted to single-group models without predictors", {
  f1 <- fit_lpa(wm, mglpa_model(1, ind, group = "group"), starts = c(1, 1), se = "none")
  f2 <- fit_lpa(wm, mglpa_model(2, ind, group = "group"), starts = c(2, 1), se = "none")
  expect_error(vlmr_test(f2, f1), "single-group")
})
