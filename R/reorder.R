#' Reorder (relabel) the profiles of a fitted model
#'
#' Latent profiles are only identified up to a permutation of their labels.
#' This function relabels the profiles of a fitted model, either with an
#' explicit permutation or automatically by the average standardized level of
#' the indicator means (`"level"`, lowest profile first) or by profile size
#' (`"size"`, largest first), and refits the model from the relabelled
#' solution so that all estimates, standard errors and posterior
#' probabilities are consistent with the new labels. The log-likelihood is
#' unchanged.
#'
#' @param fit An `mglpa_fit` object.
#' @param order Either `"level"`, `"size"`, or an integer permutation of
#'   `1:k` giving, for each new profile label, the old label it takes over.
#'   For models whose profile means are group-specific (the configural
#'   model), a list with one permutation per group may be given; the
#'   automatic rules are then applied within each group.
#' @return An `mglpa_fit` object.
#' @export
reorder_profiles <- function(fit, order = "level") {
  stopifnot(inherits(fit, "mglpa_fit"))
  model <- fit$model
  K <- model$k
  G <- fit$D$G
  if (K == 1L) return(fit)
  group_specific <- !model$equal_means
  if (is.character(order)) {
    order <- match.arg(order, c("level", "size"))
    perms <- lapply(seq_len(G), function(g) profile_order_rule(fit, order, g))
    if (!group_specific) {
      # a single ordering for all groups: average the criterion over groups
      crit <- Reduce(`+`, lapply(seq_len(G), function(g) profile_order_criterion(fit, order, g))) / G
      perms <- rep(list(order(crit)), G)
    }
  } else if (is.list(order)) {
    if (!group_specific) stop("A list of permutations is only meaningful when profile means are group-specific.")
    if (length(order) != G) stop("`order` must have one permutation per group.")
    perms <- lapply(order, check_perm, K = K)
  } else {
    perms <- rep(list(check_perm(order, K)), G)
  }
  if (all(vapply(perms, function(p) identical(p, seq_len(K)), logical(1)))) return(fit)

  new_coef <- permute_values(fit$coefficients, perms, fit)
  new_fixed <- if (length(fit$fixed)) permute_values(fit$fixed, perms, fit) else NULL
  maxiter <- if (is.null(fit$starts$maxiter)) 500 else fit$starts$maxiter
  tol <- if (is.null(fit$starts$tol)) 1e-6 else fit$starts$tol
  out <- fit_lpa_core(fit$D, model, starts = c(0, 0), stiterations = 1, maxiter = maxiter,
                      tol = tol, fixed = new_fixed, start_values = new_coef,
                      se = fit$se_type, mplus_convention = fit$mplus_convention,
                      data_frame = fit$data_frame)
  if (abs(out$loglik - fit$loglik) > 1e-3)
    warning("The log-likelihood changed after relabelling (", round(fit$loglik, 3), " -> ",
            round(out$loglik, 3), ").", call. = FALSE)
  # The refit is a single optimisation started from the relabelled solution,
  # so its own replication record is meaningless: keep the random-start record
  # of the original fit and rebuild the estimation warnings from it (the
  # convergence and information-matrix status come from the refit).
  out$starts <- fit$starts
  out$time <- fit$time
  out$warnings <- estimation_warnings(
    converged = out$converged, iterations = out$iterations, tol = tol, K = K,
    n_final = length(fit$starts$final_ll), reps = fit$starts$best_ll_replications,
    n_failed = if (is.null(fit$starts$n_failed)) 0L else fit$starts$n_failed,
    singular = isTRUE(out$information_singular))
  out
}

check_perm <- function(p, K) {
  p <- as.integer(p)
  if (length(p) != K || !setequal(p, seq_len(K))) stop("`order` must be a permutation of 1:", K)
  p
}

profile_order_criterion <- function(fit, rule, g) {
  K <- fit$model$k
  if (rule == "size") return(-fit$class_counts$by_group_model[g, ])
  pm <- profile_means(fit, standardized = TRUE)
  pm <- pm[pm$group == g & pm$variable %in% fit$model$indicators, ]
  vapply(seq_len(K), function(k) mean(pm$z[pm$class == k]), numeric(1))
}

profile_order_rule <- function(fit, rule, g) order(profile_order_criterion(fit, rule, g))

# Apply per-group permutations to a vector of values named by parameter label.
# perms[[g]][k_new] = k_old. The structure of every label (type, variable,
# profile, group, and whether it is shared across groups or profiles) is read
# from the parameter table of the fit, so variable names are never parsed.
# Membership logits (intercepts and slopes) are re-expressed relative to the
# new reference profile: theta_new[k'] = theta_old[perm[k']] - theta_old[perm[K]].
permute_values <- function(vals, perms, fit) {
  pt <- fit$params
  K <- fit$model$k
  labs <- names(vals)
  out <- vals
  first <- pt[match(labs, pt$label), ]
  shared_g <- vapply(labs, function(l) length(unique(pt$group[pt$label == l])) > 1L, logical(1))
  shared_k <- vapply(labs, function(l) length(unique(pt$class[pt$label == l])) > 1L, logical(1))
  g_eff <- ifelse(shared_g, NA_integer_, first$group)
  key <- function(type, variable, class, group)
    paste(type, variable, class, ifelse(is.na(group), "shared", group), sep = "|")
  keys <- key(first$type, first$variable, first$class, g_eff)
  lab_by_key <- stats::setNames(labs, keys)
  is_logit <- first$type %in% c("intercept", "slope")
  for (i in which(!is_logit & !shared_k)) {
    perm <- perms[[if (is.na(g_eff[i])) 1L else g_eff[i]]]
    k_old <- perm[first$class[i]]
    src <- lab_by_key[key(first$type[i], first$variable[i], k_old, g_eff[i])]
    if (!is.na(src)) out[i] <- vals[src]
  }
  if (any(is_logit)) {
    stems <- unique(data.frame(type = first$type[is_logit], variable = first$variable[is_logit],
                               g = g_eff[is_logit], stringsAsFactors = FALSE))
    for (s in seq_len(nrow(stems))) {
      g <- stems$g[s]
      perm <- perms[[if (is.na(g)) 1L else g]]
      v_old <- numeric(K)
      for (k in seq_len(K - 1L)) {
        src <- lab_by_key[key(stems$type[s], stems$variable[s], k, g)]
        if (!is.na(src)) v_old[k] <- vals[src]
      }
      v_new <- v_old[perm] - v_old[perm[K]]
      for (k in seq_len(K - 1L)) {
        dst <- lab_by_key[key(stems$type[s], stems$variable[s], k, g)]
        if (!is.na(dst)) out[dst] <- v_new[k]
      }
    }
  }
  out
}
