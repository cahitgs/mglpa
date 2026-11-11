#' Profile means in raw and standardized metric
#'
#' Returns the estimated indicator (and outcome) means of each profile, with
#' standard errors, and their standardized counterparts
#' `z = (mean - sample mean) / sample SD`, where the sample statistics are
#' computed on the data used in the analysis (pooled across groups by
#' default, as in the profile plots of the MG-LPA Analyzer).
#'
#' @param fit An `mglpa_fit` object.
#' @param standardized Logical; add standardized means.
#' @param standardize_by `"pooled"` uses the sample statistics of the whole
#'   analysed sample; `"group"` standardizes within each known group.
#' @return A data frame with columns `group`, `group_level`, `class`,
#'   `variable`, `mean`, `se`, `ci_lower`, `ci_upper`, `variance`, `z` and
#'   `proportion` (model-based profile proportion in the group).
#' @export
profile_means <- function(fit, standardized = TRUE, standardize_by = c("pooled", "group")) {
  stopifnot(inherits(fit, "mglpa_fit"))
  standardize_by <- match.arg(standardize_by)
  pt <- fit$params
  m <- pt[pt$type == "mean", ]
  v <- pt[pt$type == "variance", ]
  out <- data.frame(group = m$group, group_level = as.character(fit$data$group_levels[m$group]),
                    class = m$class, variable = m$variable, mean = m$est, se = m$se,
                    ci_lower = m$ci_lower, ci_upper = m$ci_upper,
                    variance = v$est[match(paste(m$variable, m$class, m$group),
                                           paste(v$variable, v$class, v$group))],
                    stringsAsFactors = FALSE)
  if (standardized) {
    if (standardize_by == "pooled") {
      ss <- fit$sample_stats$pooled
      idx <- match(out$variable, ss$variable)
      out$z <- (out$mean - ss$mean[idx]) / ss$sd[idx]
    } else {
      ss <- fit$sample_stats$by_group
      idx <- match(paste(out$group, out$variable), paste(ss$group, ss$variable))
      out$z <- (out$mean - ss$mean[idx]) / ss$sd[idx]
    }
  }
  pr <- fit$class_counts$by_group_model
  out$proportion <- pr[cbind(out$group, out$class)]
  rownames(out) <- NULL
  out
}

#' Profile proportions by group
#'
#' @param fit An `mglpa_fit` object.
#' @param type `"model"` (model-implied), `"posterior"` (average posterior
#'   probabilities) or `"most_likely"` (modal assignment).
#' @return A matrix with one row per group and one column per profile.
#' @export
profile_proportions <- function(fit, type = c("model", "posterior", "most_likely")) {
  type <- match.arg(type)
  switch(type, model = fit$class_counts$by_group_model,
         posterior = fit$class_counts$by_group_posterior,
         most_likely = fit$class_counts$by_group_most_likely)
}
