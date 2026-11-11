#' APA-style tables of the results
#'
#' Formatting helpers that turn the result objects of the package into
#' publication-ready tables (data frames of character columns).
#'
#' @param x A result object (`mglpa_enumeration`, `mglpa_similarity`,
#'   `mglpa_predictive`, `mglpa_explanatory` or `mglpa_classification`).
#' @param digits Number of decimals.
#' @return A data frame.
#' @name apa_tables
NULL

fmt <- function(x, digits = 3) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "< .001", sub("^0", "", formatC(p, format = "f", digits = 3))))

#' @rdname apa_tables
#' @export
apa_table <- function(x, digits = 3) UseMethod("apa_table")

#' @rdname apa_tables
#' @export
apa_table.mglpa_enumeration <- function(x, digits = 3) {
  tab <- x$table
  out <- data.frame(
    Group = if (identical(unique(tab$group), "all")) NULL else tab$group,
    Model = paste0(tab$k, "-profile"),
    LL = fmt(tab$LL, digits), fp = tab$npar, AIC = fmt(tab$AIC, digits),
    BIC = fmt(tab$BIC, digits), CAIC = fmt(tab$CAIC, digits), ABIC = fmt(tab$ABIC, digits),
    Entropy = fmt(tab$entropy, digits),
    check.names = FALSE, stringsAsFactors = FALSE)
  tests <- character()
  if ("aLMR_p" %in% names(tab)) {
    out[["aLMR p"]] <- fmt_p(tab$aLMR_p)
    tests <- c(tests, "aLMR p = p-value of the Lo-Mendell-Rubin adjusted likelihood ratio test of k versus k - 1 profiles")
  }
  if ("BLRT_p" %in% names(tab)) {
    out[["BLRT p"]] <- fmt_p(tab$BLRT_p)
    tests <- c(tests, "BLRT p = p-value of the bootstrap likelihood ratio test of k versus k - 1 profiles")
  }
  out[["Smallest profile (%)"]] <- fmt(100 * tab$min_proportion, 1)
  out[["Best LL replicated"]] <- paste0(tab$replications, "/", tab$n_final)
  attr(out, "note") <- paste0(
    "LL = log-likelihood; fp = number of free parameters; AIC = Akaike information criterion; ",
    "BIC = Bayesian information criterion; CAIC = consistent AIC; ABIC = sample-size adjusted BIC",
    if (length(tests)) paste0("; ", paste(tests, collapse = "; ")) else "", ".")
  out
}

#' @rdname apa_tables
#' @export
apa_table.mglpa_similarity <- function(x, digits = 3) {
  tab <- x$table
  out <- data.frame(
    Model = paste0(toupper(substr(tab$model, 1, 1)), substr(tab$model, 2, nchar(tab$model)), " similarity"),
    LL = fmt(tab$LL, digits), fp = tab$npar, AIC = fmt(tab$AIC, digits),
    BIC = fmt(tab$BIC, digits), CAIC = fmt(tab$CAIC, digits), ABIC = fmt(tab$ABIC, digits),
    Entropy = fmt(tab$entropy, digits),
    `Delta CAIC` = fmt(tab$dCAIC, digits), `Delta BIC` = fmt(tab$dBIC, digits),
    `Delta ABIC` = fmt(tab$dABIC, digits),
    Decision = ifelse(is.na(tab$decision), "", tab$decision),
    check.names = FALSE, stringsAsFactors = FALSE)
  attr(out, "note") <- paste(
    "LL = log-likelihood; fp = free parameters; AIC = Akaike information criterion;",
    "BIC = Bayesian information criterion; CAIC = consistent AIC; ABIC = sample-size adjusted BIC.",
    "Differences are relative to the preceding model; similarity is supported when at least two of",
    "CAIC, BIC and ABIC decrease.")
  out
}

#' @rdname apa_tables
#' @export
apa_table.mglpa_predictive <- function(x, digits = 3) {
  rows <- lapply(names(x$effects), function(nm) {
    e <- x$effects[[nm]]
    data.frame(Model = paste(nm, "effects"), Group = e$group_level, Predictor = e$predictor,
               Comparison = e$comparison,
               `Coef. (SE)` = sprintf("%s (%s)%s", fmt(e$coef, digits), fmt(e$se, digits),
                                      ifelse(e$p < .05, "*", "")),
               p = fmt_p(e$p), OR = fmt(e$OR, digits),
               `95% CI of OR` = sprintf("[%s; %s]", fmt(e$OR_ci_lower, digits), fmt(e$OR_ci_upper, digits)),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "note") <- "* p < .05. Coef. = multinomial logistic regression coefficient; SE = standard error; OR = odds ratio."
  out
}

#' @rdname apa_tables
#' @export
apa_table.mglpa_explanatory <- function(x, digits = 3) {
  rows <- lapply(names(x$contrasts), function(nm) {
    e <- x$contrasts[[nm]]
    data.frame(Model = nm, Outcome = e$outcome, Type = e$type, Group = e$group_level,
               Comparison = e$comparison,
               `Difference (SE)` = sprintf("%s (%s)%s", fmt(e$estimate, digits), fmt(e$se, digits),
                                           ifelse(e$p < .05, "*", "")),
               `95% CI` = sprintf("[%s; %s]", fmt(e$ci_lower, digits), fmt(e$ci_upper, digits)),
               p = fmt_p(e$p), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "note") <- "* p < .05. Differences in outcome means estimated with the delta method."
  out
}

#' Write a table as a Word-compatible document
#'
#' Writes the table as an HTML file with a `.doc` extension that Microsoft
#' Word and LibreOffice open directly, with APA-style formatting (Times New
#' Roman, horizontal rules only). No additional software is required.
#'
#' @param table A data frame, e.g. from [apa_table()].
#' @param file Output file name (`.doc`).
#' @param title Table title (printed in italics above the table).
#' @param note Table note; defaults to the `note` attribute of `table`.
#' @param landscape Logical; landscape page orientation.
#' @return The file name, invisibly.
#' @export
write_table_doc <- function(table, file, title = NULL, note = attr(table, "note"),
                            landscape = TRUE) {
  tables <- list(table)
  names(tables) <- if (is.null(title)) "" else title
  write_tables_doc(tables, file, landscape = landscape, notes = list(note))
}

# Several tables in one Word-compatible document. `tables` is a named list of
# data frames (names are the titles); notes default to the "note" attributes.
write_tables_doc <- function(tables, file, landscape = TRUE, notes = NULL) {
  esc <- function(s) {
    s <- gsub("&", "&amp;", s, fixed = TRUE)
    s <- gsub("<", "&lt;", s, fixed = TRUE)
    gsub(">", "&gt;", s, fixed = TRUE)
  }
  one_table <- function(table, title, note) {
    head_cells <- paste0("<th>", esc(names(table)), "</th>", collapse = "")
    body <- apply(table, 1, function(r) paste0("<tr>", paste0("<td", ifelse(seq_along(r) == 1,
      ' style="text-align:left"', ""), ">", esc(as.character(r)), "</td>", collapse = ""), "</tr>"))
    c(if (!is.null(title) && nzchar(title)) paste0("<p class='title'>", esc(title), "</p>"),
      "<table><thead><tr>", head_cells, "</tr></thead><tbody>", body, "</tbody></table>",
      if (!is.null(note)) paste0("<p class='note'><i>Note.</i> ", esc(note), "</p>"))
  }
  blocks <- unlist(lapply(seq_along(tables), function(i) {
    note <- if (!is.null(notes) && length(notes) >= i) notes[[i]] else attr(tables[[i]], "note")
    c(one_table(tables[[i]], names(tables)[i], note),
      if (i < length(tables)) "<br clear=all style='page-break-before:always'>")
  }))
  html <- c(
    "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:w='urn:schemas-microsoft-com:office:word' xmlns='http://www.w3.org/TR/REC-html40'>",
    "<head><meta charset='utf-8'><title>Tables</title>",
    "<style>",
    sprintf("@page Section1 { size: %s; margin: 2.54cm; %s }",
            if (landscape) "29.7cm 21cm" else "21cm 29.7cm",
            if (landscape) "mso-page-orientation: landscape;" else ""),
    "div.Section1 { page: Section1; }",
    "body { font-family: 'Times New Roman', serif; font-size: 12pt; }",
    "p.title { font-style: italic; margin-bottom: 8pt; }",
    "table { border-collapse: collapse; width: 100%; }",
    "th { border-top: 1.5pt solid black; border-bottom: 1pt solid black; padding: 4pt 6pt; font-weight: normal; text-align: center; }",
    "td { padding: 3pt 6pt; text-align: center; }",
    "tr:last-child td { border-bottom: 1.5pt solid black; }",
    "p.note { font-size: 11pt; margin-top: 8pt; }",
    "</style></head><body><div class='Section1'>",
    blocks,
    "</div></body></html>")
  con <- file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(c(intToUtf8(0xFEFF), html), con)
  invisible(file)
}
