# Additional Mplus reference runs produced with Mplus 9.1.1 on the example data
# (see dev/mplus_reference_extra.R): class-varying variances, three known
# groups and missing indicator values handled by FIML.
test_that("class-varying variances, three groups and missing data reproduce Mplus", {
  skip_on_cran()
  extra <- mplus_ref[mplus_ref$step == "extra", ]
  skip_if(nrow(extra) == 0, "no extra Mplus reference values")
  ref <- function(model, column) extra[[column]][extra$model == model]
  wm3 <- wm
  wm3$group3 <- cut(wm3$education, c(-Inf, 14, 16, Inf), labels = FALSE)
  wm_miss <- read.csv(system.file("extdata", "mglpa_wm_missing.csv", package = "mglpa"))
  check <- function(name, data, model, starts = c(200, 20)) {
    f <- fit_lpa(data, model, starts = starts, stiterations = 20, seed = 7, se = "none")
    expect_equal(f$loglik, ref(name, "ll"), tolerance = 5e-6, label = paste(name, "LL"))
    expect_equal(f$npar, ref(name, "npar"), label = paste(name, "npar"))
    expect_equal(f$fit[["BIC"]], ref(name, "bic"), tolerance = 5e-6, label = paste(name, "BIC"))
    # Mplus 8 reports the pattern-based entropy, Mplus 9 the profile-specific one
    if (!is.na(ref(name, "entropy")))
      expect_equal(f$entropy, ref(name, "entropy"), tolerance = 3e-3, label = paste(name, "entropy"))
    if (!is.na(ref(name, "entropy_c")))
      expect_equal(f$entropy_profiles, ref(name, "entropy_c"), tolerance = 3e-3,
                   label = paste(name, "entropy (C)"))
    invisible(f)
  }
  check("m1_g1_k3_classvarying", wm1, mglpa_model(3, ind, variances = "class_varying"))
  check("m2_configural_classvarying", wm,
        mglpa_model(3, ind, group = "group", variances = "class_varying"))
  check("m3_dispersion_classvarying", wm,
        mglpa_model(3, ind, group = "group", variances = "class_varying", equal_means = TRUE,
                    equal_variances = TRUE))
  f4 <- check("m4_structural_3groups", wm3, mglpa_model(3, ind, group = "group3", equal_means = TRUE))
  expect_equal(f4$D$G, 3)
  expect_equal(f4$npar, f4$npar_free + 2)   # two known-class proportions counted by Mplus
  check("m5_configural_3groups", wm3, mglpa_model(3, ind, group = "group3"))
  f6 <- check("m6_g1_k3_missing", wm_miss[wm_miss$group == 1, ], mglpa_model(3, ind))
  expect_true(any(is.na(wm_miss[wm_miss$group == 1, ind])))
  expect_equal(f6$data$N, 500)
})
