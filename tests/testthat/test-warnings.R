test_that("non-convergence is detected, extended and reported", {
  # tol = 0 can never be met, so the EM stops at the iteration limit
  expect_warning(
    fit <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(3, 1), stiterations = 1, maxiter = 1,
                   tol = 0, seed = 1, se = "none"),
    "DID NOT TERMINATE NORMALLY")
  expect_false(fit$converged)
  expect_true(any(grepl("DID NOT TERMINATE", fit$warnings)))
  # the best solution received additional iterations (4 x maxiter) before giving up
  expect_equal(fit$starts$extra_iterations, 4)
  expect_equal(fit$iterations, 5)
  expect_output(print(fit), "Estimation warnings")
  # with a reachable criterion the same model converges and carries no warning
  fit2 <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(3, 1), seed = 1, se = "none")
  expect_true(fit2$converged)
  expect_length(fit2$warnings, 0)
})

test_that("a replication warning is given when the best log-likelihood is not replicated", {
  # a single final-stage optimisation cannot be judged, no warning
  expect_silent(fit_lpa(wm1, mglpa_model(2, ind), starts = c(5, 1), seed = 2, se = "none"))
  # force two different local maxima by running only one EM iteration per start
  res <- withCallingHandlers(
    fit_lpa(wm1, mglpa_model(3, ind), starts = c(6, 6), stiterations = 1, maxiter = 1,
            seed = 3, se = "none"),
    warning = function(w) invokeRestart("muffleWarning"))
  expect_true(res$starts$best_ll_replications < 6)
  expect_true(length(res$warnings) >= 1)
})

test_that("relabelling the profiles keeps the replication record and its warnings", {
  fit <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(10, 3), seed = 4, se = "none")
  expect_equal(fit$starts$maxiter, 500)
  expect_equal(fit$starts$tol, 1e-6)
  expect_false(isTRUE(fit$information_singular))
  # pretend that the best log-likelihood was found in only one of the three
  # final-stage optimisations
  fit$starts$best_ll_replications <- 1L
  fit$warnings <- mglpa:::estimation_warnings(TRUE, fit$iterations, 1e-6, 3L, 3L, 1L)
  expect_true(any(grepl("NOT REPLICATED", fit$warnings)))
  r <- reorder_profiles(fit, c(2, 3, 1))   # a refit from the relabelled solution
  expect_equal(r$loglik, fit$loglik, tolerance = 1e-6)
  expect_equal(r$starts$best_ll_replications, 1L)
  expect_true(any(grepl("NOT REPLICATED", r$warnings)))
  expect_true(r$converged)
  # the weaker note is kept as well, and full replication gives no note
  fit$starts$best_ll_replications <- 2L
  expect_true(any(grepl("replicated only 2 times", reorder_profiles(fit, c(2, 3, 1))$warnings)))
  fit$starts$best_ll_replications <- 3L
  expect_length(reorder_profiles(fit, c(2, 3, 1))$warnings, 0)
  # the enumeration and similarity tables take their notes from the relabelled fits
  expect_match(mglpa:::fit_row(r)$notes, "NOT REPLICATED")
})

test_that("tables and apps carry the estimation notes", {
  enum <- enumerate_profiles(wm1, ind, k = 1:2, starts = c(4, 2), seed = 1, verbose = FALSE)
  expect_true("notes" %in% names(enum$table))
  expect_true(all(enum$table$converged))
  expect_equal(enum$table$notes, c("", ""))
})
