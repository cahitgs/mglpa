test_that("classification tables reproduce the Mplus output", {
  skip_on_cran()
  f4 <- fit_lpa(wm1, mglpa_model(4, ind), starts = c(500, 50), seed = 44, se = "none",
                verbose = FALSE)
  expect_equal(f4$loglik, ref_value("enumeration", "ll", group = 1, k = 4), tolerance = 1e-6)
  cq <- classification_quality(f4)
  expect_s3_class(cq, "mglpa_classification")
  expect_equal(cq$entropy, 0.933, tolerance = 1e-3)
  # Mplus orders the profiles differently; match them by the most likely counts
  perm <- match(c(80, 23, 219, 178), cq$counts[[1]])
  expect_false(anyNA(perm))
  A <- unname(cq$average_posterior[[1]][perm, perm])
  mplus_A <- matrix(c(0.881, 0.111, 0.008, 0.000,
                      0.223, 0.773, 0.004, 0.000,
                      0.008, 0.002, 0.988, 0.003,
                      0.000, 0.000, 0.007, 0.993), 4, byrow = TRUE)
  expect_lt(max(abs(A - mplus_A)), 1.5e-3)
  expect_equal(unname(rowSums(A)), rep(1, 4))
  P <- unname(cq$classification_probabilities[[1]][perm, perm])
  mplus_P <- matrix(c(0.912, 0.066, 0.021, 0.000,
                      0.329, 0.658, 0.014, 0.000,
                      0.003, 0.000, 0.991, 0.006,
                      0.000, 0.000, 0.003, 0.997), 4, byrow = TRUE)
  expect_lt(max(abs(P - mplus_P)), 1.5e-3)
  expect_equal(unname(rowSums(P)), rep(1, 4))
  expect_output(print(cq), "Average posterior")
  tab <- apa_table(cq)
  expect_equal(nrow(tab), 4)
  expect_true(all(c("Most likely profile", "n", "Profile 4") %in% names(tab)))
})

test_that("classification quality is computed within groups", {
  fit <- fit_lpa(wm, mglpa_model(3, ind, group = "group", equal_means = TRUE), starts = c(10, 2),
                 seed = 5, se = "none")
  cq <- classification_quality(fit)
  expect_equal(names(cq$average_posterior), c("1", "2"))
  expect_equal(names(cq$entropy_by_group), c("1", "2"))
  expect_equal(unname(rowSums(cq$average_posterior[["2"]])), rep(1, 3))
  expect_equal(sum(cq$counts[["1"]]), 500)
  tab <- apa_table(cq)
  expect_equal(names(tab)[1], "Group")
  expect_equal(nrow(tab), 6)
  expect_output(print(cq), "Group 2")
  f1 <- fit_lpa(wm1, mglpa_model(1, ind), starts = c(1, 1), se = "none")
  expect_output(print(classification_quality(f1)), "not defined")
})
