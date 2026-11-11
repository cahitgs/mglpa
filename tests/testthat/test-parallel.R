test_that("arguments given as variables reach the parallel worker processes", {
  skip_on_cran()
  # The functions that run on the workers are created inside fit_lpa() and
  # blrt(); arguments that are still unevaluated when they are sent would be
  # looked up in the caller's environment, and variables of the global
  # environment do not exist on the workers.
  assign("mglpa_test_starts", c(3, 2), envir = globalenv())
  assign("mglpa_test_iter", 5, envir = globalenv())
  assign("mglpa_test_tol", 1e-5, envir = globalenv())
  on.exit(rm(list = c("mglpa_test_starts", "mglpa_test_iter", "mglpa_test_tol"),
             envir = globalenv()), add = TRUE)
  fit <- suppressWarnings(eval(bquote(
    mglpa::fit_lpa(.(wm1), mglpa::mglpa_model(3, .(ind)), starts = mglpa_test_starts,
                   stiterations = mglpa_test_iter, maxiter = 200 + mglpa_test_iter,
                   tol = mglpa_test_tol, seed = 2, se = "none", cores = 2)), globalenv()))
  expect_s3_class(fit, "mglpa_fit")
  expect_equal(length(fit$starts$final_ll), 2)
  expect_equal(fit$starts$maxiter, 205)
  expect_true(is.finite(fit$loglik))

  f2 <- fit_lpa(wm1, mglpa_model(2, ind), starts = c(10, 2), seed = 1, se = "none")
  res <- suppressWarnings(eval(bquote(
    mglpa::blrt(.(fit), .(f2), B = 4, starts = mglpa_test_starts, seed = 5, verbose = FALSE,
                cores = 2)), globalenv()))
  expect_equal(res$successful, 4)
  expect_true(all(is.finite(res$bootstrap)))
  expect_equal(res$draws$n_final, rep(2, 4))
  expect_error(blrt(fit, f2, B = 2, starts = 5), "two non-negative integers")
})
