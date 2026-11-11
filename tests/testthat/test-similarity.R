test_that("the similarity sequence reproduces the Mplus results of the article", {
  skip_on_cran()
  sim <- similarity_tests(wm, ind, group = "group", k = 3, starts = c(60, 10), seed = 21,
                          verbose = FALSE)
  for (nm in c("configural", "structural", "dispersion", "distributional")) {
    expect_equal(sim$fits[[nm]]$loglik, ref_value(nm, "ll"), tolerance = 1e-5, label = nm)
    expect_equal(sim$fits[[nm]]$npar, ref_value(nm, "npar"), label = nm)
    expect_equal(sim$fits[[nm]]$fit[["BIC"]], ref_value(nm, "bic"), tolerance = 1e-5, label = nm)
    expect_equal(sim$fits[[nm]]$entropy, 0.981, tolerance = 2e-3, label = nm)
  }
  expect_equal(sim$table$decision, c(NA, "supported", "supported", "rejected"))
  expect_equal(sim$retained, "dispersion")
  # profile proportions by group of the dispersion model (Mplus: 20.9/43.7/35.4 and 44.7/40.6/14.7)
  pr <- profile_proportions(sim$fits$dispersion)
  expect_equal(unname(sort(pr[1, ])), c(0.209, 0.354, 0.437), tolerance = 5e-3)
  expect_equal(unname(sort(pr[2, ])), c(0.147, 0.406, 0.447), tolerance = 5e-3)
  # profiles are ordered by level: profile 1 lowest
  pm <- profile_means(sim$fits$dispersion)
  lev <- vapply(1:3, function(k) mean(pm$z[pm$class == k & pm$group == 1]), numeric(1))
  expect_true(all(diff(lev) > 0))
  expect_s3_class(apa_table(sim), "data.frame")
  expect_output(print(sim), "Most constrained model supported: dispersion")

  # ---- predictive similarity -------------------------------------------------
  pred <- predictive_similarity(sim, predictors = c("education", "cog_reserve"), seed = 22,
                                verbose = FALSE)
  expect_equal(pred$fits$free$loglik, ref_value("predictive", "ll", model = "predictive freely freevar"),
               tolerance = 1e-5)
  expect_equal(pred$fits$equal$loglik, ref_value("predictive", "ll", model = "predictive equal freevar"),
               tolerance = 1e-5)
  expect_equal(pred$fits$free$npar, 18)
  expect_equal(pred$fits$equal$npar, 14)
  expect_equal(pred$decision$decision, "supported")
  eff <- pred$effects$free
  expect_equal(nrow(eff), 2 * 2 * 3)   # groups x predictors x profile pairs
  expect_true(all(eff$OR > 0))
  expect_s3_class(apa_table(pred), "data.frame")

  # ---- explanatory similarity ------------------------------------------------
  expl <- explanatory_similarity(sim, outcomes = c("task_perf", "everyday_func"), seed = 23,
                                 verbose = FALSE)
  ref_names <- c(free_means_free_var = "explanatory freely freevar",
                 equal_means_free_var = "explanatory equal freevar",
                 free_means_equal_var = "explanatory freely fixedvar",
                 equal_means_equal_var = "explanatory equal fixedvar")
  for (nm in names(ref_names)) {
    expect_equal(expl$fits[[nm]]$loglik, ref_value("explanatory", "ll", model = ref_names[[nm]]),
                 tolerance = 1e-5, label = nm)
    expect_equal(expl$fits[[nm]]$npar, ref_value("explanatory", "npar", model = ref_names[[nm]]),
                 label = nm)
  }
  expect_equal(expl$decisions$equal$decision, "rejected")
  expect_equal(expl$decisions$free$decision, "rejected")
  # Mplus (equal means, equal variances): task_perf profile 1 - profile 2 = -13.633 (SE 0.767)
  ct <- expl$contrasts$equal_means_equal_var
  d12 <- ct[ct$outcome == "task_perf" & ct$comparison == "Profile 1 - Profile 2", ]
  expect_equal(d12$estimate, -13.633, tolerance = 1e-3)
  expect_equal(d12$se, 0.767, tolerance = 0.02)
  expect_equal(d12$ci_lower, -15.137, tolerance = 2e-3)
  expect_s3_class(apa_table(expl), "data.frame")
  expect_output(print(expl), "Explanatory similarity")
})

test_that("compare_fit applies the two-out-of-three rule", {
  mk <- function(caic, bic, abic) structure(list(fit = c(CAIC = caic, BIC = bic, ABIC = abic)))
  expect_equal(compare_fit(mk(1, 1, 1), mk(2, 2, 2))$decision, "supported")
  expect_equal(compare_fit(mk(1, 1, 3), mk(2, 2, 2))$decision, "supported")
  expect_equal(compare_fit(mk(1, 3, 3), mk(2, 2, 2))$decision, "mixed")
  expect_equal(compare_fit(mk(3, 3, 3), mk(2, 2, 2))$decision, "rejected")
})
