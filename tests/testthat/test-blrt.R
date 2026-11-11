test_that("data can be simulated from a fitted model and the BLRT runs", {
  f2 <- fit_lpa(wm1, mglpa_model(2, ind), starts = c(10, 2), seed = 1, se = "none")
  f3 <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(10, 2), seed = 1, se = "none")
  sim <- simulate_lpa(f2, n = 200)
  expect_equal(dim(sim), c(200, 6))
  expect_true(all(sim$profile %in% 1:2))
  res <- suppressWarnings(blrt(f3, f2, B = 4, starts = c(5, 2), seed = 2, verbose = FALSE))
  expect_s3_class(res, "mglpa_blrt")
  expect_equal(length(res$bootstrap), 4)
  expect_equal(res$successful, 4)
  expect_true(is.data.frame(res$draws))
  expect_true(all(c("statistic", "converged_k", "converged_km1", "replications", "n_final") %in%
                    names(res$draws)))
  expect_equal(res$draws$n_final, rep(2, 4))
  expect_true(is.logical(res$draws$converged_k))
  expect_true(res$p_upper >= res$p && res$p >= res$p_lower)
  expect_true(is.character(res$warnings))
  expect_gt(res$statistic, 0)
  expect_true(res$p >= 0 && res$p <= 1)
  expect_equal(res$p, mean(res$bootstrap >= res$statistic))
  expect_output(print(res), "Bootstrap LRT")
  res2 <- suppressWarnings(blrt(f3, f2, B = 4, starts = c(3, 1), seed = 5, verbose = FALSE, cores = 2))
  expect_equal(length(res2$bootstrap), 4)
  expect_error(blrt(f2, f3), "one more profile")
  expect_error(simulate_lpa(fit_lpa(wm, mglpa_model(2, ind, group = "group"), starts = c(2, 1),
                                    se = "none")), "single-group")
})

test_that("BLRT diagnostics flag draws that may be local maxima, as Mplus does", {
  draws <- data.frame(statistic = c(1, 2, 20, 3, NA, 25, 30),
                      converged_k = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE),
                      converged_km1 = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
                      replications = c(3, 3, 1, 1, 1, 1, 3), n_final = rep(4, 7))
  s <- mglpa:::blrt_summary(draws, stat_obs = 10)
  expect_equal(s$successful, 6)
  expect_equal(s$p, 3 / 6)
  # understated: draw 2 (k model not converged) and draw 4 (not replicated)
  # are below the observed value; draw 3 is not replicated but exceeds it
  expect_equal(s$n_flagged, 2)
  expect_equal(s$p_upper, 5 / 6)
  # overstated: draw 6 exceeds the observed value with a non-converged
  # (k - 1)-profile model
  expect_equal(s$n_flagged_high, 1)
  expect_equal(s$p_lower, 2 / 6)
  # draw 7 (k model not converged but above the observed value) is only
  # reported in the convergence note
  expect_equal(s$n_not_converged, 3)
  expect_true(any(grepl("LOCAL MAXIMA", s$warnings)))
  expect_true(any(grepl("may be overstated", s$warnings)))
  expect_true(any(grepl("3 of 6 bootstrap draws did not meet the convergence criterion", s$warnings)))
  expect_true(any(grepl("1 of 7 bootstrap draws failed", s$warnings)))
  # a non-converged draw never stays silent, even when it exceeds the observed value
  s_nc <- mglpa:::blrt_summary(draws[c(1, 7), ], 10)
  expect_equal(s_nc$n_flagged, 0)
  expect_equal(s_nc$n_not_converged, 1)
  expect_true(any(grepl("convergence criterion", s_nc$warnings)))
  s2 <- mglpa:::blrt_summary(draws[c(1, 3), ], 10)
  expect_equal(s2$n_flagged, 0)
  expect_length(s2$warnings, 0)
  # a single final-stage optimisation cannot be judged for replication
  s3 <- mglpa:::blrt_summary(data.frame(statistic = c(1, 2), converged_k = TRUE,
                                        converged_km1 = TRUE, replications = 1, n_final = 1), 10)
  expect_equal(s3$n_flagged, 0)
  expect_equal(s3$p, 0)
  # a negative statistic is a local maximum even when it was replicated
  s4 <- mglpa:::blrt_summary(data.frame(statistic = c(-5, 12), converged_k = TRUE,
                                        converged_km1 = TRUE, replications = 4, n_final = 4), 10)
  expect_equal(s4$n_flagged, 1)
  expect_equal(s4$p_upper, 1)
})

test_that("BLRT p-values are added to an enumeration table", {
  enum <- enumerate_profiles(wm1, ind, k = 1:2, starts = c(5, 2), seed = 3, verbose = FALSE)
  enum2 <- enumeration_blrt(enum, B = 3, starts = c(3, 1), seed = 4, verbose = FALSE)
  expect_true("BLRT_p" %in% names(enum2$table))
  expect_true(is.na(enum2$table$BLRT_p[enum2$table$k == 1]))
  expect_false(is.na(enum2$table$BLRT_p[enum2$table$k == 2]))
  expect_equal(which(names(enum2$table) == "BLRT_p"), which(names(enum2$table) == "aLMR_p") + 1L)
  expect_s3_class(enum2$blrt[["all 2"]], "mglpa_blrt")
  expect_output(print(enum2), "BLRT_p")
})
