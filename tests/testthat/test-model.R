test_that("model specification validates its inputs", {
  expect_s3_class(mglpa_model(2, c("a", "b")), "mglpa_model")
  expect_error(mglpa_model(0, "a"), "positive integer")
  expect_error(mglpa_model(2, character()), "non-empty")
  expect_error(mglpa_model(2, c("a", "a")), "only once")
  expect_error(mglpa_model(2, "a", equal_means = TRUE), "require a `group`")
  expect_output(print(mglpa_model(2, c("a", "b"), group = "g", predictors = "x", outcomes = "y")),
                "predictors: x")
})

test_that("parameter tables reproduce the Mplus free-parameter counts", {
  n_free <- function(...) {
    pt <- mglpa:::build_param_table(mglpa_model(3, ind, group = "group", ...), G = 2)
    length(unique(pt$label))
  }
  # Mplus adds one parameter for the known-class proportion (two groups)
  expect_equal(n_free() + 1, 45)                                             # configural
  expect_equal(n_free(equal_means = TRUE) + 1, 30)                           # structural
  expect_equal(n_free(equal_means = TRUE, equal_variances = TRUE) + 1, 25)   # dispersion
  expect_equal(n_free(equal_means = TRUE, equal_variances = TRUE,
                      equal_proportions = TRUE) + 1, 23)                     # distributional
  # class-varying variances
  expect_equal(n_free(variances = "class_varying"), 30 + 30 + 4)
})

test_that("parameter rows are laid out in array order", {
  m <- mglpa_model(2, c("a", "b", "c"), group = "g", predictors = c("x1", "x2"), outcomes = "y")
  pt <- mglpa:::build_param_table(m, G = 2)
  means <- pt[pt$type == "mean", ]
  expect_equal(means$pos, seq_len(nrow(means)))
  expect_equal(means$variable[1:4], c("a", "b", "c", "y"))
  expect_equal(means$class[1:8], rep(1:2, each = 4))
  slopes <- pt[pt$type == "slope", ]
  expect_equal(slopes$pos, seq_len(nrow(slopes)))
  expect_equal(unique(slopes$label), c("b_x1_c1_g1", "b_x2_c1_g1", "b_x1_c1_g2", "b_x2_c1_g2"))
  # outcome variances: one per outcome by default
  expect_equal(unique(pt$label[pt$type == "variance" & pt$variable == "y"]), "v_y")
})

test_that("fixed parameters are removed from the free set", {
  m <- mglpa_model(2, c("a", "b"), group = "g")
  pt <- mglpa:::build_param_table(m, G = 2)
  pt2 <- mglpa:::apply_fixed(pt, c(m_a_c1_g1 = 0.5, v_a_g1 = 2))
  expect_false(any(pt2$free[pt2$label %in% c("m_a_c1_g1", "v_a_g1")]))
  expect_equal(unique(pt2$value[pt2$label == "v_a_g1"]), 2)
  expect_error(mglpa:::apply_fixed(pt, c(nonsense = 1)), "Unknown parameter labels")
})
