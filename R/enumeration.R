#' Class enumeration: profile solutions with 1 to K profiles in each group
#'
#' Fits latent profile models with an increasing number of profiles,
#' separately within each known group (step 0 of the multi-group similarity
#' workflow), and collects the fit indices used to decide on the number of
#' profiles. Consecutive solutions are compared with the Vuong-Lo-Mendell-Rubin
#' likelihood ratio tests ([vlmr_test()], Mplus `TECH11`); the bootstrap
#' likelihood ratio test can be added with [enumeration_blrt()].
#'
#' @inheritParams mglpa_model
#' @inheritParams fit_lpa
#' @param k Integer vector of the numbers of profiles to estimate.
#' @param order_profiles How profiles are labelled after estimation: by
#'   average standardized level of the indicator means (`"level"`, lowest
#'   first), by size (`"size"`, largest first) or left as estimated
#'   (`"none"`). See [reorder_profiles()].
#' @param lrt Logical; compute the VLMR and adjusted LMR tests for each
#'   solution against the solution with one profile less.
#' @param verbose Logical; report progress.
#' @return An object of class `mglpa_enumeration` with elements `fits` (a
#'   list by group of lists by `k` of `mglpa_fit` objects), `table` (one row
#'   per model with log-likelihood, free parameters, AIC, BIC, CAIC, ABIC,
#'   entropy, the likelihood ratio statistic and the VLMR and aLMR p-values,
#'   the smallest profile proportion, and the number of final-stage
#'   optimisations that replicated the best log-likelihood) and `vlmr` (the
#'   test objects).
#' @examples
#' \donttest{
#' data(mglpa_wm)
#' enum <- enumerate_profiles(mglpa_wm,
#'   indicators = c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed"),
#'   group = "group", k = 1:4, starts = c(20, 4), seed = 1)
#' enum
#' }
#' @export
enumerate_profiles <- function(data, indicators, group = NULL, k = 1:8,
                               variances = c("class_invariant", "class_varying"),
                               starts = c(100, 10), stiterations = 10, maxiter = 500,
                               tol = 1e-6, se = c("robust", "ml", "none"),
                               order_profiles = c("level", "size", "none"), lrt = TRUE,
                               group_levels = NULL, seed = NULL, verbose = TRUE,
                               cores = 1L) {
  variances <- match.arg(variances)
  se <- match.arg(se)
  order_profiles <- match.arg(order_profiles)
  if (!is.null(seed)) set.seed(seed)
  k <- sort(unique(as.integer(k)))
  if (!is.null(group)) {
    if (is.null(data[[group]])) stop("Group variable `", group, "` not found in `data`.")
    gv <- data[[group]]
    lev <- if (is.null(group_levels)) sort(unique(gv[!is.na(gv)])) else group_levels
    subsets <- lapply(lev, function(l) data[!is.na(gv) & gv == l, , drop = FALSE])
    names(subsets) <- as.character(lev)
  } else {
    lev <- NULL
    subsets <- list(all = data)
  }
  fits <- list()
  rows <- list()
  for (gname in names(subsets)) {
    fits[[gname]] <- list()
    for (kk in k) {
      model <- mglpa_model(kk, indicators, variances = variances)
      f <- fit_lpa(subsets[[gname]], model, starts = starts, stiterations = stiterations,
                   maxiter = maxiter, tol = tol, se = se, verbose = FALSE, cores = cores)
      if (order_profiles != "none" && kk > 1L) f <- reorder_profiles(f, order_profiles)
      fits[[gname]][[as.character(kk)]] <- f
      rows[[length(rows) + 1L]] <- fit_row(f, group = gname, k = kk)
      if (verbose)
        message(sprintf("%s k = %d: LL = %.3f, BIC = %.3f, entropy = %s, best LL replicated %d/%d",
                        if (gname == "all") "" else paste0("Group ", gname), kk, f$loglik,
                        f$fit[["BIC"]], if (is.na(f$entropy)) "NA" else sprintf("%.3f", f$entropy),
                        f$starts$best_ll_replications, length(f$starts$final_ll)))
    }
  }
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  out <- structure(list(fits = fits, table = table, indicators = indicators, group = group,
                        group_levels = lev, variances = variances, k = k, data = data,
                        call = match.call()),
                   class = "mglpa_enumeration")
  if (lrt && length(k) > 1L)
    out <- tryCatch(enumeration_vlmr(out), error = function(e) {
      warning("The VLMR tests could not be computed: ", conditionMessage(e), call. = FALSE)
      out
    })
  out
}

fit_row <- function(f, group = NA, k = f$model$k, model = NA_character_) {
  props <- f$class_counts$most_likely$proportion
  data.frame(model = model, group = group, k = k, LL = f$loglik, npar = f$npar,
             AIC = f$fit[["AIC"]], BIC = f$fit[["BIC"]], CAIC = f$fit[["CAIC"]],
             ABIC = f$fit[["ABIC"]], entropy = f$entropy, min_proportion = min(props),
             replications = f$starts$best_ll_replications,
             n_final = length(f$starts$final_ll), converged = f$converged,
             notes = paste(f$warnings, collapse = " "),
             stringsAsFactors = FALSE)
}

#' @export
print.mglpa_enumeration <- function(x, digits = 3, ...) {
  cat("<mglpa_enumeration> ", length(x$k), " solutions (k = ", paste(range(x$k), collapse = "-"),
      ") in ", length(x$fits), " group(s); ", gsub("_", "-", x$variances), " variances\n", sep = "")
  tab <- x$table
  for (g in unique(tab$group)) {
    tg <- tab[tab$group == g, ]
    if (g != "all") cat("\nGroup ", g, "\n", sep = "")
    out <- data.frame(k = tg$k, LL = round(tg$LL, digits), npar = tg$npar,
                      AIC = round(tg$AIC, digits), BIC = round(tg$BIC, digits),
                      CAIC = round(tg$CAIC, digits), ABIC = round(tg$ABIC, digits),
                      entropy = round(tg$entropy, digits))
    for (col in c("VLMR_p", "aLMR_p", "BLRT_p"))
      if (col %in% names(tg)) out[[col]] <- round(tg[[col]], 3)
    out$min_prop <- round(tg$min_proportion, 3)
    out$reps <- paste0(tg$replications, "/", tg$n_final)
    flag <- function(v) ifelse(v == min(v, na.rm = TRUE), "*", "")
    out$best <- paste0(flag(tg$BIC), ifelse(tg$CAIC == min(tg$CAIC), "C", ""),
                       ifelse(tg$ABIC == min(tg$ABIC), "A", ""))
    out$conv <- ifelse(tg$converged, "yes", "NO")
    print(out, row.names = FALSE)
  }
  cat("\n* lowest BIC, C lowest CAIC, A lowest ABIC within group.",
      if ("VLMR_p" %in% names(tab))
        "VLMR_p / aLMR_p: Vuong-Lo-Mendell-Rubin tests of k versus k - 1 profiles (Mplus TECH11)."
      else "",
      if ("BLRT_p" %in% names(tab)) "BLRT_p: bootstrap likelihood ratio test (TECH14)." else "",
      "\n  Inspect the elbow plot and the substantive meaning of the profiles before retaining k.\n")
  fits <- unlist(x$fits, recursive = FALSE)
  names(fits) <- unlist(lapply(names(x$fits), function(g)
    paste0(if (g == "all") "" else paste0("Group ", g, ", "), "k = ", names(x$fits[[g]]))))
  print_estimation_notes(fits)
  invisible(x)
}

#' Fit-index table of a class enumeration
#'
#' @param x An `mglpa_enumeration` object.
#' @return A data frame.
#' @export
enumeration_table <- function(x) {
  stopifnot(inherits(x, "mglpa_enumeration"))
  x$table
}
