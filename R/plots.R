#' @importFrom rlang .data
NULL

#' Colour palettes for profile plots
#'
#' The nine palettes of the MG-LPA Analyzer web application.
#' @format A named list of character vectors of hexadecimal colours.
#' @export
mglpa_palettes <- list(
  vibrant = c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6", "#1abc9c", "#e67e22", "#34495e"),
  pastel = c("#a6cee3", "#fdbf6f", "#b2df8a", "#fb9a99", "#cab2d6", "#ffff99", "#ff7f00", "#e31a1c"),
  ocean = c("#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600"),
  warm = c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a"),
  cool = c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61"),
  nature = c("#1b7837", "#5aae61", "#a6dba0", "#d9f0d3", "#e7d4e8", "#c2a5cf"),
  colorblind = c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"),
  grayscale = c("#000000", "#404040", "#707070", "#A0A0A0", "#C0C0C0", "#E0E0E0"),
  apa = c("#0077C8", "#7F7F7F", "#000000", "#5A9BD5", "#ED7D31", "#A5A5A5")
)

palette_colours <- function(palette, n) {
  cols <- if (is.character(palette) && length(palette) == 1L && palette %in% names(mglpa_palettes))
    mglpa_palettes[[palette]] else palette
  rep_len(cols, n)
}

#' Elbow plot of the information criteria of a class enumeration
#'
#' @param x An `mglpa_enumeration` object.
#' @param indices Which criteria to plot.
#' @return A `ggplot` object.
#' @export
elbow_plot <- function(x, indices = c("AIC", "BIC", "CAIC", "ABIC")) {
  stopifnot(inherits(x, "mglpa_enumeration"))
  tab <- x$table
  long <- do.call(rbind, lapply(indices, function(ic)
    data.frame(group = tab$group, k = tab$k, index = ic, value = tab[[ic]])))
  long$index <- factor(long$index, levels = indices)
  long$group <- factor(long$group, levels = unique(tab$group),
                       labels = if (identical(unique(tab$group), "all")) "All" else
                         paste("Group", unique(tab$group)))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$k, y = .data$value, colour = .data$index,
                                          shape = .data$index, group = .data$index)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = sort(unique(long$k))) +
    ggplot2::labs(x = "Number of profiles", y = "Information criterion", colour = NULL,
                  shape = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  if (nlevels(long$group) > 1L) p <- p + ggplot2::facet_wrap(~group, scales = "free_y")
  p
}

#' Profile plot
#'
#' Plots the profile-specific means of the indicators, standardized with the
#' sample mean and standard deviation of each indicator (z scores) or in
#' their original metric. In multi-group models one panel per group is drawn
#' when the plotted means differ between groups (the configural model); when
#' they are constrained to equality (structural, dispersion and
#' distributional models) a single panel is drawn and the group-specific
#' profile proportions are reported in the labels.
#'
#' @param fit An `mglpa_fit` object.
#' @param type `"bar"` (profiles on the x axis, one bar per indicator),
#'   `"line"` (indicators on the x axis, one line per profile) or `"radar"`.
#' @param standardized Logical; plot z scores (default) or raw means.
#' @param variables Variables to plot; defaults to the profile indicators.
#' @param palette Name of a palette in [mglpa_palettes] or a vector of colours.
#' @param ylim Optional limits of the value axis (radial axis for `"radar"`).
#' @param title Optional plot title.
#' @param show_proportions Logical; report the model-based profile
#'   proportions in the labels.
#' @param standardize_by Passed to [profile_means()].
#' @param facet `"auto"` (one panel per group only when the plotted means
#'   differ between groups), `"group"` (always one panel per group) or
#'   `"none"` (a single panel showing the means of the first group).
#' @return A `ggplot` object.
#' @export
profile_plot <- function(fit, type = c("bar", "line", "radar"), standardized = TRUE,
                         variables = NULL, palette = "vibrant", ylim = NULL, title = NULL,
                         show_proportions = TRUE, standardize_by = "pooled",
                         facet = c("auto", "group", "none")) {
  stopifnot(inherits(fit, "mglpa_fit"))
  type <- match.arg(type)
  facet <- match.arg(facet)
  if (is.null(variables)) variables <- fit$model$indicators
  pm <- profile_means(fit, standardized = TRUE, standardize_by = standardize_by)
  pm <- pm[pm$variable %in% variables, ]
  pm$value <- if (standardized) pm$z else pm$mean
  pm$variable <- factor(pm$variable, levels = variables)
  G <- length(fit$data$group_levels)
  K <- fit$model$k
  lev <- as.character(fit$data$group_levels)
  pr <- fit$class_counts$by_group_model   # G x K model-based proportions

  # Do the plotted means coincide across groups (equality constraints)?
  key <- paste(pm$class, pm$variable)
  v1 <- pm$value[pm$group == 1L][order(key[pm$group == 1L])]
  same_means <- G > 1L && all(vapply(seq_len(G)[-1L], function(g) {
    vg <- pm$value[pm$group == g][order(key[pm$group == g])]
    length(vg) == length(v1) && max(abs(vg - v1)) < 1e-8
  }, logical(1)))
  single <- facet == "none" || (facet == "auto" && (G == 1L || same_means))
  if (single && G > 1L) pm <- pm[pm$group == 1L, ]
  same_props <- G > 1L && max(abs(sweep(pr, 2, pr[1L, ]))) < 1e-6
  glab <- if (all(grepl("^[0-9]+$", lev))) paste0("G", lev) else lev
  pct <- function(x) sprintf("%.1f%%", 100 * x)
  # proportion text for profile k: one value, or one value per group
  prop_text <- function(k) {
    if (G == 1L || same_props) return(pct(pr[1L, k]))
    vapply(k, function(kk) paste(paste(glab, pct(pr[, kk])), collapse = " | "), "")
  }
  group_subtitle <- function() paste(vapply(seq_len(G), function(g)
    sprintf("Group %s: %s", lev[g], paste(sprintf("P%d %s", seq_len(K), pct(pr[g, ])), collapse = ", ")),
    ""), collapse = "   |   ")
  subtitle <- NULL
  pm$group_label <- factor(paste("Group", pm$group_level), levels = paste("Group", lev))

  if (type == "bar") {
    if (!single) {
      # one panel per group: labels carry the proportion of that group
      pm$profile <- if (show_proportions) sprintf("Profile %d\n(%s)", pm$class, pct(pm$proportion)) else
        sprintf("Profile %d", pm$class)
      pm$profile <- factor(pm$profile, levels = unique(pm$profile[order(pm$group, pm$class)]))
    } else if (show_proportions && (G <= 3L || same_props)) {
      pm$profile <- sprintf("Profile %d\n(%s)", pm$class, prop_text(pm$class))
      pm$profile <- factor(pm$profile, levels = unique(pm$profile[order(pm$class)]))
    } else {
      pm$profile <- factor(sprintf("Profile %d", pm$class),
                           levels = sprintf("Profile %d", seq_len(K)))
      if (show_proportions) subtitle <- group_subtitle()
    }
  } else {
    if (show_proportions && (G == 1L || same_props)) {
      pm$profile <- sprintf("Profile %d (%s)", pm$class, pct(pr[1L, pm$class]))
    } else {
      pm$profile <- sprintf("Profile %d", pm$class)
      if (show_proportions) subtitle <- group_subtitle()
    }
    pm$profile <- factor(pm$profile, levels = unique(pm$profile[order(pm$class)]))
  }

  ylab <- if (standardized) "Standardized score (z)" else "Mean"
  if (type == "bar") {
    p <- ggplot2::ggplot(pm, ggplot2::aes(x = .data$profile, y = .data$value, fill = .data$variable)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      ggplot2::scale_fill_manual(values = palette_colours(palette, length(variables))) +
      ggplot2::labs(x = NULL, y = ylab, fill = NULL)
  } else if (type == "line") {
    p <- ggplot2::ggplot(pm, ggplot2::aes(x = .data$variable, y = .data$value,
                                          colour = .data$profile, group = .data$profile)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_colour_manual(values = palette_colours(palette, nlevels(pm$profile))) +
      ggplot2::labs(x = NULL, y = ylab, colour = NULL)
  } else {
    # Radar: the indicators are placed at equal angles on a numeric axis whose
    # end point (p + 1) coincides with the first indicator, which closes the
    # polygons without showing a duplicated label.
    n_var <- length(variables)
    pm$xnum <- as.integer(pm$variable)
    closing <- pm[pm$xnum == 1L, ]
    closing$xnum <- n_var + 1L
    pm_r <- rbind(pm, closing)
    r_lim <- if (!is.null(ylim)) ylim else
      c(floor((min(pm$value) - 0.5) * 2) / 2, ceiling((max(pm$value) + 0.5) * 2) / 2)
    p <- ggplot2::ggplot(pm_r, ggplot2::aes(x = .data$xnum, y = .data$value,
                                            colour = .data$profile, group = .data$profile)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_x_continuous(breaks = seq_len(n_var), labels = variables,
                                  limits = c(1, n_var + 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = r_lim) +
      ggplot2::scale_colour_manual(values = palette_colours(palette, nlevels(pm$profile))) +
      ggplot2::coord_radial(start = 0, expand = FALSE, inner.radius = 0.05) +
      ggplot2::labs(x = NULL, y = NULL, colour = NULL)
  }
  if (!is.null(ylim) && type != "radar") p <- p + ggplot2::coord_cartesian(ylim = ylim)
  if (!single) p <- p + ggplot2::facet_wrap(~group_label, scales = if (type == "bar") "free_x" else "fixed")
  p + ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
                   axis.text.x = ggplot2::element_text(angle = if (type == "bar") 0 else 20,
                                                       hjust = if (type == "bar") 0.5 else 1))
}

#' @export
plot.mglpa_fit <- function(x, ...) profile_plot(x, ...)

#' @export
plot.mglpa_enumeration <- function(x, ...) elbow_plot(x, ...)
