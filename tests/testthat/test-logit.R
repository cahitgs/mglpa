test_that("the multinomial logit solver recovers known coefficients", {
  set.seed(1)
  N <- 4000
  x <- rnorm(N)
  g <- rep(1:2, each = N / 2)
  Xd <- cbind(as.numeric(g == 1), as.numeric(g == 2), x * (g == 1), x * (g == 2))
  th_true <- matrix(c(0.5, -0.3, 1.0, -0.8,
                      -0.2, 0.4, 0.3, 0.9), ncol = 2)
  eta <- Xd %*% th_true
  pi <- exp(cbind(eta, 0)); pi <- pi / rowSums(pi)
  cls <- apply(pi, 1, function(p) sample(3, 1, prob = p))
  tau <- matrix(0, N, 3); tau[cbind(seq_len(N), cls)] <- 1
  th0 <- matrix(0, 4, 2)
  mask <- matrix(TRUE, 4, 2)
  th <- mglpa:::logit_newton(tau, Xd, th0, mask)
  expect_equal(th, th_true, tolerance = 0.15)
  # gradient at the solution is zero
  logpi <- mglpa:::softmax_ref(Xd %*% th)
  grad <- crossprod(Xd, tau[, -3] - exp(logpi)[, -3])
  expect_lt(max(abs(grad)), 1e-5)
  # fixed entries stay fixed
  mask2 <- mask; mask2[1, 1] <- FALSE
  th1 <- th0; th1[1, 1] <- 0.123
  th2 <- mglpa:::logit_newton(tau, Xd, th1, mask2)
  expect_equal(th2[1, 1], 0.123)
})

test_that("predictors of profile membership enter the model", {
  fit <- fit_lpa(wm, mglpa_model(2, ind, group = "group", predictors = "education",
                                  predictor_effects = "equal"), starts = c(10, 2), seed = 8)
  expect_true("b_education_c1" %in% names(coef(fit)))
  # Mplus convention: known-class proportion + saturated normal for the predictor
  expect_equal(fit$npar, fit$npar_free + 1 + 2)
  expect_equal(fit$loglik, fit$loglik_mixture + fit$loglik_known_groups + fit$loglik_predictors)
  eff <- predictor_effects(fit)
  expect_equal(nrow(eff), 1)
  expect_equal(eff$OR, exp(eff$coef))
  fit0 <- fit_lpa(wm, mglpa_model(2, ind, group = "group", predictors = "education",
                                   predictor_effects = "equal"), starts = c(10, 2), seed = 8,
                  mplus_convention = FALSE)
  expect_equal(fit0$npar, fit0$npar_free)
})
