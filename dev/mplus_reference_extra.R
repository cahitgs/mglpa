# Generate additional Mplus reference runs (class-varying variances, three
# groups, missing data) with the local Mplus installation and compare them
# with the mglpa engine. Writes mglpa/dev/mplus_extra/*.inp/.out and prints a
# comparison table; the parsed reference values are appended to
# inst/extdata/mplus_reference_fits.csv by the second part of the script.
args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[1] else "run"
pkgload::load_all("mglpa", compile = FALSE, quiet = TRUE)
mplus <- "C:/Program Files/Mplus/Mplus.exe"
dir <- "mglpa/dev/mplus_extra"
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
d <- read.csv("mglpa/inst/extdata/mglpa_wm_dataset_n1000.csv")
ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
mplus_names <- c(id = "id", group = "group", group3 = "group3", verbal_wm = "verbal_m",
                 visual_wm = "visual_m", updating = "updating", inhibition = "inhibitn",
                 proc_speed = "proc_spd", education = "educatin", cog_reserve = "cog_rese",
                 task_perf = "task_pef", everyday_func = "everydac")
mind <- unname(mplus_names[ind])

# Three groups: years of education cut at 14 and 16 ---------------------------
d$group3 <- cut(d$education, c(-Inf, 14, 16, Inf), labels = FALSE)
# Missing data: 8% of every indicator missing completely at random -------------
set.seed(2026)
dm <- d
for (v in ind) dm[[v]][stats::runif(nrow(dm)) < 0.08] <- NA
allna <- rowSums(is.na(dm[, ind])) == length(ind)
dm[allna, ind[1]] <- d[allna, ind[1]]
write.csv(dm[, c("id", "group", "group3", ind, "education", "cog_reserve", "task_perf", "everyday_func")],
          "mglpa/inst/extdata/mglpa_wm_missing.csv", row.names = FALSE)

write_dat <- function(df, file) {
  x <- df[, names(mplus_names)]
  x[is.na(x)] <- 999
  write.table(x, file, row.names = FALSE, col.names = FALSE, sep = "\t")
}
write_dat(d, file.path(dir, "wm3.dat"))
write_dat(dm, file.path(dir, "wm_miss.dat"))

# Mplus truncates input lines longer than 90 characters
wrap_names <- function(x, width = 70) {
  out <- character(); line <- ""
  for (v in x) {
    if (nchar(line) + nchar(v) + 1 > width) { out <- c(out, line); line <- v } else
      line <- if (nzchar(line)) paste(line, v) else v
  }
  paste(c(out, line), collapse = "\n    ")
}
names_line <- wrap_names(mplus_names)
vars_line <- paste(mind, collapse = " ")
lab <- function(prefix, idx) paste0(prefix, idx)
analysis <- "ANALYSIS:\n  TYPE = MIXTURE;\n  ESTIMATOR = MLR;\n  STARTS = 1000 100;\n  STITERATIONS = 20;\n  PROCESSORS = 8;\n"
block_lines <- function(means_lab, var_lab) {
  sprintf("  [%s] (%s-%s);\n  %s (%s-%s);", vars_line, lab(means_lab, 1), lab(means_lab, 5),
          vars_line, lab(var_lab, 1), lab(var_lab, 5))
}
mg_model <- function(G, K, means, variances) {
  # means: "group" (m{g}_{p}_) or "equal" (m{p}_); variances: "gk", "k", "g", "equal"
  logit <- if (G > 1) paste(sprintf("  c#%d on %s;", seq_len(K - 1),
                                    paste(sprintf("cg#%d", seq_len(G - 1)), collapse = " ")),
                            collapse = "\n") else ""
  blocks <- character()
  for (g in seq_len(G)) for (p in seq_len(K)) {
    ml <- if (means == "group") sprintf("m%d_%d_", g, p) else sprintf("m%d_", p)
    vl <- switch(variances, gk = sprintf("v%d_%d_", g, p), k = sprintf("v%d_", p),
                 g = sprintf("vg%d_", g), equal = "v_")
    blocks <- c(blocks, sprintf("  %%cg#%d.c#%d%%\n%s", g, p, block_lines(ml, vl)))
  }
  paste(c("MODEL:\n  %OVERALL%", logit, blocks), collapse = "\n")
}
inputs <- list(
  m1_g1_k3_classvarying = list(
    data = "wm3.dat", useobs = "group EQ 1", known = NULL, K = 3,
    model = paste(c("MODEL:\n  %OVERALL%",
                    sprintf("  %%c#%d%%\n  [%s];\n  %s;", 1:3, vars_line, vars_line)), collapse = "\n")),
  m2_configural_classvarying = list(data = "wm3.dat", known = "cg (group = 1 group = 2)", G = 2, K = 3,
                                    model = mg_model(2, 3, "group", "gk")),
  m3_dispersion_classvarying = list(data = "wm3.dat", known = "cg (group = 1 group = 2)", G = 2, K = 3,
                                    model = mg_model(2, 3, "equal", "k")),
  m4_structural_3groups = list(data = "wm3.dat", known = "cg (group3 = 1 group3 = 2 group3 = 3)",
                               G = 3, K = 3, model = mg_model(3, 3, "equal", "g")),
  m5_configural_3groups = list(data = "wm3.dat", known = "cg (group3 = 1 group3 = 2 group3 = 3)",
                               G = 3, K = 3, model = mg_model(3, 3, "group", "g")),
  m6_g1_k3_missing = list(
    data = "wm_miss.dat", useobs = "group EQ 1", missing = TRUE, K = 3,
    model = paste(c("MODEL:\n  %OVERALL%", sprintf("  %%c#%d%%\n  [%s];", 1:3, vars_line)),
                  collapse = "\n"))
)
for (nm in names(inputs)) {
  s <- inputs[[nm]]
  var <- c(sprintf("  NAMES = %s;", names_line), sprintf("  USEVARIABLES = %s;", vars_line),
           if (!is.null(s$useobs)) sprintf("  USEOBS = %s;", s$useobs),
           if (isTRUE(s$missing)) "  MISSING = ALL (999);",
           if (!is.null(s$known)) sprintf("  KNOWNCLASS = %s;", s$known),
           if (!is.null(s$known)) sprintf("  CLASSES = cg (%d) c (%d);", s$G, s$K) else
             sprintf("  CLASSES = c (%d);", s$K))
  inp <- c(sprintf("TITLE: %s", nm), sprintf("DATA: FILE IS %s;", s$data), "VARIABLE:", var, "",
           analysis, s$model, "", "OUTPUT: TECH1;")
  writeLines(inp, file.path(dir, paste0(nm, ".inp")))
}

parse_mplus <- function(file) {
  s <- paste(readLines(file, warn = FALSE), collapse = "\n")
  get <- function(p) { m <- regmatches(s, regexec(p, s))[[1]]; if (length(m) < 2) NA_real_ else as.numeric(m[2]) }
  # Mplus 8 prints one pattern-based entropy under CLASSIFICATION QUALITY;
  # Mplus 9 prints C-specific results (entropy for CG and for C)
  list(ll = get("H0 Value\\s+(-?[0-9.]+)"), npar = get("Number of Free Parameters\\s+([0-9]+)"),
       aic = get("Akaike \\(AIC\\)\\s+([0-9.]+)"), bic = get("Bayesian \\(BIC\\)\\s+([0-9.]+)"),
       abic = get("Sample-Size Adjusted BIC\\s+([0-9.]+)"),
       entropy = if (grepl("Classification Quality for C\\b", s))
         get("CLASSIFICATION QUALITY\\s+Entropy\\s+([0-9.]+)") else get("Entropy\\s+([0-9.]+)"),
       entropy_c = get("Classification Quality for C\\s+Entropy\\s+([0-9.]+)"),
       n = get("Number of observations\\s+([0-9]+)"),
       normal = grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", s),
       replicated = grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", s))
}

if (what == "run") {
  for (nm in names(inputs)) {
    t0 <- Sys.time()
    withr::with_dir(dir, system2(mplus, shQuote(paste0(nm, ".inp")), stdout = FALSE, stderr = FALSE))
    r <- parse_mplus(file.path(dir, paste0(nm, ".out")))
    cat(sprintf("%-28s Mplus LL = %12.3f npar = %3d BIC = %11.3f entropy = %.3f N = %d normal = %s replicated = %s  (%.0fs)\n",
                nm, r$ll, r$npar, r$bic, r$entropy, r$n, r$normal, r$replicated,
                as.numeric(Sys.time() - t0, units = "secs")))
    flush.console()
  }
}

if (what == "compare") {
  d1 <- d[d$group == 1, ]
  models <- list(
    m1_g1_k3_classvarying = list(data = d1, model = mglpa_model(3, ind, variances = "class_varying")),
    m2_configural_classvarying = list(data = d, model = mglpa_model(3, ind, group = "group", variances = "class_varying")),
    m3_dispersion_classvarying = list(data = d, model = mglpa_model(3, ind, group = "group", variances = "class_varying",
                                                                    equal_means = TRUE, equal_variances = TRUE)),
    m4_structural_3groups = list(data = d, model = mglpa_model(3, ind, group = "group3", equal_means = TRUE)),
    m5_configural_3groups = list(data = d, model = mglpa_model(3, ind, group = "group3")),
    m6_g1_k3_missing = list(data = dm[dm$group == 1, ], model = mglpa_model(3, ind))
  )
  rows <- list()
  for (nm in names(models)) {
    r <- parse_mplus(file.path(dir, paste0(nm, ".out")))
    t0 <- Sys.time()
    f <- fit_lpa(models[[nm]]$data, models[[nm]]$model, starts = c(1000, 100), stiterations = 20,
                 seed = 42, se = "none", cores = 4)
    cat(sprintf("%-28s Mplus LL %12.3f | mglpa LL %12.3f | diff %+8.4f | npar %d/%d | entropy(patterns) %s/%.3f | entropy(C) %s/%.3f | N %d/%d | reps %d/%d | %.0fs\n",
                nm, r$ll, f$loglik, f$loglik - r$ll, r$npar, f$npar,
                if (is.na(r$entropy)) "NA" else sprintf("%.3f", r$entropy), f$entropy,
                if (is.na(r$entropy_c)) "NA" else sprintf("%.3f", r$entropy_c), f$entropy_profiles,
                r$n, f$data$N, f$starts$best_ll_replications, length(f$starts$final_ll),
                as.numeric(Sys.time() - t0, units = "secs")))
    flush.console()
    rows[[nm]] <- data.frame(step = "extra", model = nm, group = "", k = 3,
                             variances = models[[nm]]$model$variances, ll = r$ll, npar = r$npar,
                             aic = r$aic, bic = r$bic, caic = r$bic + r$npar, abic = r$abic,
                             entropy = r$entropy, almr_p = NA, blrt_p = NA, proportions = "",
                             best_ll_replications = NA, entropy_c = r$entropy_c,
                             mplus_version = "9.1.1", stringsAsFactors = FALSE)
  }
  ref <- read.csv("mglpa/inst/extdata/mplus_reference_fits.csv", stringsAsFactors = FALSE)
  ref <- ref[ref$step != "extra", ]
  if (!"entropy_c" %in% names(ref)) ref$entropy_c <- NA_real_
  if (!"mplus_version" %in% names(ref)) ref$mplus_version <- "8.3"
  ref <- rbind(ref, do.call(rbind, rows))
  write.csv(ref, "mglpa/inst/extdata/mplus_reference_fits.csv", row.names = FALSE, na = "")
  cat("reference table updated:", nrow(ref), "rows\n")
}
