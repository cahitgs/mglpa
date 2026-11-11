# Run the same EM iterations from the same starting values with both engines.
em_both <- function(model, data, iters = 25, seed = 1, fixed = NULL) {
  D <- mglpa:::prepare_data(data, model)
  pt <- mglpa:::apply_fixed(mglpa:::build_param_table(model, D$G), fixed)
  info_r <- mglpa:::make_info(pt, D, model, engine = "r")
  info_c <- mglpa:::make_info(pt, D, model, engine = "auto")
  set.seed(seed)
  par <- mglpa:::vals_to_par(mglpa:::random_start(info_r$pt, D, info_r, model), info_r)
  list(r = mglpa:::em_iterate(par, D, info_r, model, maxiter = iters, tol = 1e-6, check = FALSE),
       cpp = mglpa:::em_iterate(par, D, info_c, model, maxiter = iters, tol = 1e-6, check = FALSE),
       use_cpp = info_c$use_cpp)
}

cases <- list(
  single = list(mglpa_model(3, ind), wm1),
  varying = list(mglpa_model(2, ind, variances = "class_varying"), wm1),
  configural = list(mglpa_model(3, ind, group = "group"), wm),
  structural = list(mglpa_model(3, ind, group = "group", equal_means = TRUE), wm),
  dispersion = list(mglpa_model(3, ind, group = "group", equal_means = TRUE,
                                equal_variances = TRUE), wm),
  distributional = list(mglpa_model(3, ind, group = "group", equal_means = TRUE,
                                    equal_variances = TRUE, equal_proportions = TRUE), wm),
  one = list(mglpa_model(1, ind, group = "group"), wm)
)

test_that("the compiled EM iterations reproduce the pure-R iterations exactly", {
  for (nm in names(cases)) {
    b <- em_both(cases[[nm]][[1]], cases[[nm]][[2]])
    expect_true(b$use_cpp, label = nm)
    expect_equal(b$cpp$ll, b$r$ll, tolerance = 1e-10, label = nm)
    expect_equal(b$cpp$par$mu, b$r$par$mu, tolerance = 1e-9, label = nm)
    expect_equal(b$cpp$par$sig2, b$r$par$sig2, tolerance = 1e-9, label = nm)
    if (!is.null(b$r$par$alpha))
      expect_equal(b$cpp$par$alpha, b$r$par$alpha, tolerance = 1e-9, label = nm)
  }
})

test_that("full fits agree between engines", {
  for (nm in c("single", "dispersion", "one")) {
    f_r <- fit_lpa(cases[[nm]][[2]], cases[[nm]][[1]], starts = c(5, 2), seed = 3, engine = "r",
                   se = "none")
    f_c <- fit_lpa(cases[[nm]][[2]], cases[[nm]][[1]], starts = c(5, 2), seed = 3, se = "none")
    expect_true(f_c$info$use_cpp, label = nm)
    expect_false(f_r$info$use_cpp, label = nm)
    expect_equal(f_c$loglik, f_r$loglik, tolerance = 1e-7, label = nm)
    expect_equal(f_c$coefficients, f_r$coefficients, tolerance = 1e-4, label = nm)
    expect_lt(max(abs(f_c$posterior - f_r$posterior)), 1e-4, label = nm)
    expect_equal(dim(f_c$class_counts$by_group_model), c(f_c$D$G, f_c$model$k))
  }
})

test_that("the compiled EM handles missing data and fixed parameters", {
  set.seed(7)
  w <- wm1
  w$visual_wm[sample(nrow(w), 60)] <- NA
  w$proc_speed[sample(nrow(w), 40)] <- NA
  b <- em_both(mglpa_model(3, ind), w)
  expect_equal(b$cpp$ll, b$r$ll, tolerance = 1e-10)
  expect_equal(b$cpp$par$mu, b$r$par$mu, tolerance = 1e-9)
  base <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(10, 2), seed = 3, se = "none")
  fx <- fixed_from_fit(base, "mean")
  b2 <- em_both(mglpa_model(3, ind), wm1, fixed = fx)
  expect_true(b2$use_cpp)
  expect_equal(b2$cpp$ll, b2$r$ll, tolerance = 1e-10)
  expect_equal(as.numeric(b2$cpp$par$mu), unname(fx[paste0("m_", rep(ind, 3), "_c", rep(1:3, each = 5), "_g1")]))
  f_c <- fit_lpa(wm1, mglpa_model(3, ind), starts = c(3, 1), fixed = fx, seed = 4, se = "none")
  expect_equal(f_c$npar, 7)
})

test_that("models with predictors or fixed logits use the R engine", {
  f <- fit_lpa(wm, mglpa_model(2, ind, group = "group", predictors = "education"),
               starts = c(2, 1), seed = 1, se = "none")
  expect_false(f$info$use_cpp)
  base <- fit_lpa(wm, mglpa_model(2, ind, group = "group"), starts = c(2, 1), seed = 1, se = "none")
  f2 <- fit_lpa(wm, mglpa_model(2, ind, group = "group"), starts = c(1, 1),
                fixed = fixed_from_fit(base, "intercept"), se = "none")
  expect_false(f2$info$use_cpp)
  old <- options(mglpa.engine = "r")
  on.exit(options(old))
  f3 <- fit_lpa(wm1, mglpa_model(2, ind), starts = c(1, 1), seed = 1, se = "none")
  expect_false(f3$info$use_cpp)
})
