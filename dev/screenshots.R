# Screenshots of the Shiny application for the README (man/figures/app-*.png).
# The application is started from the installed package in a background R
# process and driven in headless Chrome with the chromote package, so the
# script also serves as an end-to-end check of the user interface.
# Usage (from the directory that contains the package folder):
#   Rscript mglpa/dev/screenshots.R          # settings used for the README
#   Rscript mglpa/dev/screenshots.R quick    # few random starts, for a trial run
args <- commandArgs(trailingOnly = TRUE)
quick <- length(args) > 0 && args[1] == "quick"
out_dir <- file.path("mglpa", "man", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
port <- 3939L
url <- sprintf("http://127.0.0.1:%d", port)
cores <- max(1L, min(8L, parallel::detectCores() - 1L))
# The class enumeration of the README uses the defaults of the application
# (1 to 6 profiles, 100 initial and 10 final random starts, seed 1): two
# solutions then carry an estimation warning, which shows the notes and the
# replication badges at work.
set <- if (quick) {
  list(kmax = 4, starts = c(50, 10), sim_starts = c(50, 10))
} else {
  list(kmax = 6, starts = c(100, 10), sim_starts = c(500, 50))
}

# ---- Start the application ----------------------------------------------------
app <- callr::r_bg(function(port) {
  mglpa::mglpa_app(port = port, host = "127.0.0.1", launch.browser = FALSE)
}, args = list(port = port))
up <- FALSE
for (i in 1:120) {
  up <- tryCatch({ suppressWarnings(readLines(url, n = 1, warn = FALSE)); TRUE },
                 error = function(e) FALSE)
  if (up) break
  if (!app$is_alive()) stop("The application stopped: ", paste(app$read_all_error(), collapse = "\n"))
  Sys.sleep(0.5)
}
if (!up) stop("The application did not start.")

# ---- Browser helpers ----------------------------------------------------------
width <- 1280
b <- chromote::ChromoteSession$new(width = width, height = 900)
js <- function(code) {
  r <- b$Runtime$evaluate(code, returnByValue = TRUE, awaitPromise = TRUE)
  if (!is.null(r$exceptionDetails))
    stop("JavaScript error: ", r$exceptionDetails$text, " ", r$exceptionDetails$exception$description)
  r$result$value
}
q <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
# Wait until the Shiny session has been idle for two seconds
wait_idle <- function(timeout = 3600, settle = 1.5) {
  t0 <- Sys.time()
  idle <- 0
  Sys.sleep(0.6)
  repeat {
    busy <- js("document.documentElement.classList.contains('shiny-busy')")
    idle <- if (isTRUE(busy)) 0 else idle + 1
    if (idle >= 4) break
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout)
      stop("Timed out waiting for the application.")
    Sys.sleep(0.5)
  }
  Sys.sleep(settle)
}
# Set an input through its Shiny input binding (updates the widget as well)
set_input <- function(id, value) {
  js(sprintf(paste0("(function(){var el=document.getElementById(%s);",
                    "var b=$(el).data('shiny-input-binding');b.setValue(el,%s);",
                    "$(el).trigger('change');return true;})()"), q(id), q(value)))
  Sys.sleep(0.4)
}
click <- function(id) js(sprintf("document.getElementById(%s).click(); true", q(id)))
open_tab <- function(title) {
  js(sprintf("document.querySelector(%s).click(); true",
             q(sprintf("a.nav-link[data-value=\"%s\"]", title))))
  wait_idle()
}
selected <- function(id) js(sprintf("JSON.stringify($(%s).val())", q(paste0("#", id))))
set_viewport <- function(w) {
  width <<- w
  b$Emulation$setDeviceMetricsOverride(width = w, height = 900, deviceScaleFactor = 1, mobile = FALSE)
  wait_idle()
}
# Collapse the settings sidebar of the open tab so that wide tables fit
collapse_sidebar <- function() {
  js(paste0("(function(){var t=document.querySelector('.tab-pane.active button.collapse-toggle');",
            "if(t&&t.getAttribute('aria-expanded')==='true')t.click();return true;})()"))
  Sys.sleep(1)
  wait_idle()
}
dismiss_notifications <- function()
  js("document.querySelectorAll('.shiny-notification-close').forEach(function(e){e.click();}); true")
# Tables of the open tab that are wider than their card (they would be cut)
cut_tables <- function() unlist(js(paste0(
  "(function(){var out=[];document.querySelectorAll('.tab-pane.active .datatables').forEach(",
  "function(w){var s=w.querySelector('.dataTables_scrollBody, .dt-scroll-body');",
  "if(s&&s.scrollWidth>s.clientWidth+1)out.push(w.id+' ('+s.scrollWidth+' > '+s.clientWidth+' px)');});",
  "return out;})()")))
# Document coordinates (left, top, right, bottom) of an element, by default
# of the card that holds it
box <- function(selector, card = TRUE) unlist(js(sprintf(paste0(
  "(function(){var el=document.querySelector(%s);var c=%s;",
  "var r=c.getBoundingClientRect();return [r.left+window.scrollX,r.top+window.scrollY,",
  "r.right+window.scrollX,r.bottom+window.scrollY];})()"), q(selector),
  if (card) "el.closest('.card')||el" else "el")))
# Screenshot from the top of the page down to the lowest of the elements in
# `until`, or of the area spanned by the cards of `cards`
shot <- function(name, until = NULL, cards = NULL, elements = NULL, pad = 14) {
  dismiss_notifications()
  Sys.sleep(0.5)
  file <- file.path(out_dir, paste0("app-", name, ".png"))
  if (!is.null(until) || !is.null(elements)) {
    bottom <- max(c(vapply(until, function(s) box(s)[4], numeric(1)),
                    vapply(elements, function(s) box(s, card = FALSE)[4], numeric(1))))
    clip <- c(0, 0, width, ceiling(bottom) + pad)
  } else {
    bx <- vapply(cards, box, numeric(4))
    clip <- c(max(0, floor(min(bx[1, ])) - pad), floor(min(bx[2, ])) - pad,
              ceiling(max(bx[3, ]) - min(bx[1, ])) + 2 * pad,
              ceiling(max(bx[4, ]) - min(bx[2, ])) + 2 * pad)
  }
  b$screenshot(file, cliprect = clip, scale = 1, delay = 0.5)
  cat(sprintf("%-40s %4.0f x %4.0f px  %4.0f KB\n", file, clip[3], clip[4], file.size(file) / 1024))
  cut <- cut_tables()
  if (length(cut)) cat("  note: tables wider than their card:", paste(cut, collapse = ", "), "\n")
}

run <- function() {
  b$Page$navigate(url)
  for (i in 1:120) {
    ready <- tryCatch(js("!!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected())"),
                      error = function(e) FALSE)
    if (isTRUE(ready)) break
    Sys.sleep(0.5)
  }
  wait_idle()

  # 1. Data: the example button loads the data and selects the variables --------
  click("use_example")
  wait_idle()
  cat("after the example button: indicators = ", selected("indicators"), ", group = ",
      selected("group"), ", predictors = ", selected("predictors"), ", outcomes = ",
      selected("outcomes"), "\n", sep = "")
  if (identical(selected("indicators"), "[]") || identical(selected("group"), "\"\""))
    stop("The example button did not select the variables.")
  set_input("cores", cores)
  shot("data", until = "#data_preview", elements = ".tab-pane.active aside.sidebar")

  # 2. Class enumeration --------------------------------------------------------
  open_tab("2. Enumeration")
  set_input("enum_k", c(1, set$kmax))
  set_input("enum_starts1", set$starts[1])
  set_input("enum_starts2", set$starts[2])
  click("run_enum")
  wait_idle()
  set_input("enum_plot_group", "1")
  set_input("enum_plot_k", "3")
  wait_idle()
  shot("profiles", cards = c("#enum_profile_plot", "#enum_classif"))
  set_viewport(1540)
  collapse_sidebar()
  cat("estimation notes of the enumeration:", js("document.querySelectorAll('#enum_notes li').length"),
      "entries\n")
  shot("enumeration", until = "#elbow")
  # The bootstrap likelihood ratio test is not shown in the README; a short
  # run checks that its p-values arrive in the table.
  set_input("blrt_B", 10)
  set_input("blrt_starts1", 10)
  set_input("blrt_starts2", 2)
  click("run_blrt")
  wait_idle()
  if (!grepl("BLRT_p", js("$('#enum_table').text()"), fixed = TRUE))
    stop("The BLRT p-values did not appear in the fit table.")
  cat("BLRT p-values shown in the fit table\n")

  # 3. Similarity tests ---------------------------------------------------------
  open_tab("3. Similarity tests")
  set_input("sim_starts1", set$sim_starts[1])
  set_input("sim_starts2", set$sim_starts[2])
  click("run_sim")
  wait_idle()
  set_input("sim_plot_type", "line")
  set_viewport(1770)
  collapse_sidebar()
  shot("similarity", until = "#sim_profile_plot")

  # 4. Predictive similarity ----------------------------------------------------
  open_tab("4. Predictive")
  click("run_pred")
  wait_idle()
  collapse_sidebar()
  shot("predictive", until = "#pred_effects_free")

  # 5. Explanatory similarity ---------------------------------------------------
  open_tab("5. Explanatory")
  set_input("expl_var", "equal")
  click("run_expl")
  wait_idle()
  collapse_sidebar()
  shot("explanatory", until = "#expl_contrasts")
}

tryCatch(run(), finally = {
  try(b$close(), silent = TRUE)
  app$kill()
})
cat("done\n")
