# Validation of the EM engine against the Mplus reference fits stored in the
# example project of the MG-LPA Analyzer (example/mglpa_project.json).
# Usage: Rscript mglpa/dev/validate_engine.R [quick|enum|sim|predexpl]
args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[1] else "quick"
pkgload::load_all("mglpa", quiet = TRUE)
d <- read.csv("mglpa/inst/extdata/mglpa_wm_dataset_n1000.csv")
ref <- read.csv("mglpa/inst/extdata/mplus_reference_fits.csv", stringsAsFactors = FALSE)
ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")

report <- function(name, fit, ll_ref, npar_ref) {
  cat(sprintf("%-28s LL = %12.3f (Mplus %12.3f, diff %+8.4f)  npar = %d (Mplus %d)  ent = %s  reps = %d/%d  conv = %s  %.1fs\n",
              name, fit$loglik, ll_ref, fit$loglik - ll_ref, fit$npar, npar_ref,
              if (is.na(fit$entropy)) "NA" else sprintf("%.3f", fit$entropy),
              fit$starts$best_ll_replications, length(fit$starts$final_ll), fit$converged, fit$time))
  flush.console()
}

if (what == "quick") {
  m <- mglpa_model(3, ind)
  f <- fit_lpa(d[d$group == 1, ], m, starts = c(50, 10), seed = 1)
  print(f)
  r <- ref[ref$step == "enumeration" & ref$group == 1 & ref$k == 3, ]
  report("Group 1 k=3", f, r$ll, r$npar)
  print(parameters(f)[1:12, ])
  cat("Mplus proportions:", r$proportions, "\n")
  print(f$class_counts$most_likely)
  print(f$class_counts$classification)
}

if (what == "enum") {
  for (g in 1:2) for (k in 1:8) {
    m <- mglpa_model(k, ind)
    f <- fit_lpa(d[d$group == g, ], m, starts = c(400, 40), stiterations = 10, seed = 100 * g + k)
    r <- ref[ref$step == "enumeration" & ref$group == g & ref$k == k, ]
    report(sprintf("Group %d k=%d", g, k), f, r$ll, r$npar)
    cat("   Mplus props:", r$proportions, " | ours:",
        paste(round(sort(f$class_counts$most_likely$proportion), 3), collapse = ";"),
        " | Mplus reps:", r$best_ll_replications, "\n")
  }
}

if (what == "sim") {
  cfg <- list(configural = c(FALSE, FALSE, FALSE), structural = c(TRUE, FALSE, FALSE),
              dispersion = c(TRUE, TRUE, FALSE), distributional = c(TRUE, TRUE, TRUE))
  for (nm in names(cfg)) {
    cc <- cfg[[nm]]
    m <- mglpa_model(3, ind, group = "group", equal_means = cc[1], equal_variances = cc[2],
                     equal_proportions = cc[3])
    f <- fit_lpa(d, m, starts = c(200, 20), seed = 7)
    r <- ref[ref$step == nm, ]
    report(nm, f, r$ll, r$npar)
    print(round(f$class_counts$by_group_model, 3))
  }
}

if (what == "predexpl") {
  m_disp <- mglpa_model(3, ind, group = "group", equal_means = TRUE, equal_variances = TRUE)
  disp <- fit_lpa(d, m_disp, starts = c(200, 20), seed = 7)
  report("dispersion", disp, ref$ll[ref$step == "dispersion"], ref$npar[ref$step == "dispersion"])
  fx <- fixed_from_fit(disp, c("mean", "variance"), digits = 5)
  X <- as.matrix(d[, c("education", "cog_reserve")]); N <- nrow(X); q <- ncol(X)
  S <- crossprod(scale(X, scale = FALSE)) / N
  ll_x <- -N / 2 * (q * log(2 * pi) + log(det(S)) + q)
  cat("Saturated MVN log-likelihood of the predictors:", ll_x, "\n")
  for (eff in c("free", "equal")) {
    m <- mglpa_model(3, ind, group = "group", equal_means = TRUE, equal_variances = TRUE,
                     predictors = c("education", "cog_reserve"), predictor_effects = eff)
    f <- fit_lpa(d, m, starts = c(20, 5), fixed = fx, start_values = disp, seed = 3)
    nm <- if (eff == "free") "predictive freely freevar" else "predictive equal freevar"
    r <- ref[ref$model == nm, ][1, ]
    cat(sprintf("%s: ours (Mplus convention, incl. MVN of predictors = %.3f) = %.3f ; Mplus = %.3f ; npar ours %d vs Mplus %d\n",
                nm, ll_x, f$loglik, r$ll, f$npar, r$npar))
    print(parameters(f, "slope")[, c("variable", "class", "group_level", "label", "est", "se", "p")])
  }
  fx2 <- fixed_from_fit(disp, c("mean", "variance", "intercept"), digits = 5)
  for (om in c("free", "equal")) for (ov in c("free", "equal")) {
    m <- mglpa_model(3, ind, group = "group", equal_means = TRUE, equal_variances = TRUE,
                     outcomes = c("task_perf", "everyday_func"), outcome_means = om,
                     outcome_variances = ov)
    f <- fit_lpa(d, m, starts = c(20, 5), fixed = fx2, seed = 3)
    nm <- sprintf("explanatory %s %s", if (om == "free") "freely" else "equal",
                  if (ov == "free") "freevar" else "fixedvar")
    r <- ref[ref$model == nm, ]
    report(nm, f, r$ll, r$npar)
    if (om == "equal" && ov == "equal")
      print(parameters(f, c("mean", "variance"))[parameters(f)$variable %in% c("task_perf", "everyday_func"), c("type", "variable", "class", "group_level", "label", "est", "se")][1:10, ])
  }
}
