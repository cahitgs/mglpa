#' Launch the MG-LPA Shiny application
#'
#' A graphical interface to the full workflow of the package: data import,
#' class enumeration, the four similarity tests, profile plots, predictive
#' and explanatory similarity, APA-style tables and project files. All
#' estimation is done in R by [fit_lpa()]; no external software is needed.
#'
#' @param data Optional data frame loaded when the application starts.
#' @param ... Passed to [shiny::runApp()], e.g. `launch.browser`.
#' @return Called for its side effect (runs the application).
#' @examples
#' if (interactive()) mglpa_app()
#' @export
mglpa_app <- function(data = NULL, ...) {
  for (p in c("shiny", "bslib", "DT")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package `", p, "` is required for the app. Install it with install.packages(\"", p, "\").",
           call. = FALSE)
  }
  app <- shiny::shinyApp(ui = app_ui(), server = app_server(data))
  shiny::runApp(app, ...)
}

# ------------------------------------------------------------------------------
# UI ---------------------------------------------------------------------------
# ------------------------------------------------------------------------------

app_ui <- function() {
  theme <- bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#1f4e79",
                           base_font = bslib::font_google("DM Sans"),
                           heading_font = bslib::font_google("Syne"))
  starts_inputs <- function(id, s1 = 100, s2 = 10) {
    shiny::tagList(
      shiny::numericInput(paste0(id, "_starts1"), "Random starts (initial stage)", s1, min = 0, step = 10),
      shiny::numericInput(paste0(id, "_starts2"), "Final-stage optimisations", s2, min = 1, step = 1),
      shiny::numericInput(paste0(id, "_stiter"), "Initial-stage iterations", 10, min = 1, step = 1),
      shiny::numericInput(paste0(id, "_seed"), "Random seed", 1, min = 1, step = 1)
    )
  }
  # Plots keep a fixed height and a comfortable maximum width instead of being
  # stretched to the full browser window.
  # The explicit width matters: inside the flex container of a card body a
  # centred box without a width shrinks to its content, the plot then reports a
  # zero width and Shiny leaks a graphics device on every failed render.
  # Cards that hold tables use a non-fillable body: in a fillable (flexbox)
  # card body a DT widget may be shrunk below the height of its table, which
  # then overflows onto the elements below it.
  table_card <- function(title, ..., full_screen = FALSE)
    bslib::card(bslib::card_header(title), bslib::card_body(fillable = FALSE, ...),
                full_screen = full_screen)
  plot_box <- function(id, height = "440px") {
    shiny::div(style = "width: 100%; max-width: 1100px; margin-left: auto; margin-right: auto;",
               shiny::plotOutput(id, height = height, width = "100%"))
  }
  bslib::page_navbar(
    title = "MG-LPA Analyzer (R)",
    theme = theme,
    id = "nav",
    fillable = FALSE,
    bslib::nav_panel(
      "1. Data",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::fileInput("file", "Upload a CSV file", accept = c(".csv", ".txt", ".dat")),
          shiny::selectInput("dec", "Decimal mark", c("." = ".", "," = ",")),
          shiny::textInput("missing_code", "Missing value code (optional)", ""),
          shiny::actionButton("use_example", "Use the example data set", class = "btn-secondary"),
          shiny::hr(),
          shiny::numericInput("cores", "CPU cores for random starts", 1, min = 1,
                              max = max(1L, parallel::detectCores(logical = TRUE)), step = 1),
          shiny::hr(),
          shiny::selectizeInput("indicators", "Profile indicators", choices = NULL, multiple = TRUE),
          shiny::selectizeInput("group", "Grouping variable (optional)", choices = NULL, multiple = FALSE),
          shiny::selectizeInput("predictors", "Predictors of profile membership (optional)",
                                choices = NULL, multiple = TRUE),
          shiny::selectizeInput("outcomes", "Distal outcomes (optional)", choices = NULL, multiple = TRUE)
        ),
        table_card("Data", shiny::uiOutput("data_summary"), DT::DTOutput("data_preview")),
        bslib::card(
          bslib::card_header("Workflow"),
          shiny::HTML(paste0(
            "<p>The application implements the multi-group tests of profile similarity of ",
            "Morin, Meyer, Creusier and Bi&eacute;try (2016): <b>class enumeration</b> in each group, ",
            "<b>configural, structural, dispersion and distributional</b> similarity, and ",
            "<b>predictive</b> and <b>explanatory</b> similarity. Models are estimated in R by ",
            "maximum likelihood (EM algorithm with random starts, robust standard errors). ",
            "Information criteria and decisions follow the same rules as the Mplus-based tool.</p>"))
        )
      )
    ),
    bslib::nav_panel(
      "2. Enumeration",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::sliderInput("enum_k", "Number of profiles", min = 1, max = 10, value = c(1, 6), step = 1),
          shiny::selectInput("enum_var", "Within-profile variances",
                             c("Equal across profiles (class-invariant)" = "class_invariant",
                               "Profile-specific (class-varying)" = "class_varying")),
          starts_inputs("enum"),
          shiny::actionButton("run_enum", "Run class enumeration", class = "btn-primary w-100"),
          shiny::hr(),
          shiny::downloadButton("dl_enum_csv", "Table (CSV)", class = "btn-sm"),
          shiny::downloadButton("dl_enum_doc", "Table (Word)", class = "btn-sm"),
          shiny::downloadButton("dl_elbow_png", "Elbow plot (PNG)", class = "btn-sm"),
          shiny::hr(),
          shiny::tags$b("Bootstrap likelihood ratio test (TECH14)"),
          shiny::numericInput("blrt_B", "Bootstrap draws", 100, min = 5, max = 5000, step = 5),
          bslib::layout_columns(
            shiny::numericInput("blrt_starts1", "k-profile model: starts", 20, min = 1, step = 5),
            shiny::numericInput("blrt_starts2", "final optimisations", 4, min = 1, step = 1)
          ),
          shiny::actionButton("run_blrt", "Run BLRT for all solutions", class = "btn-secondary w-100")
        ),
        table_card("Fit indices", DT::DTOutput("enum_table")),
        shiny::uiOutput("enum_notes"),
        bslib::card(bslib::card_header("Elbow plot"), full_screen = TRUE, plot_box("elbow", "420px")),
        bslib::card(
          bslib::card_header("Profile plot of a solution"),
          full_screen = TRUE,
          bslib::layout_columns(
            shiny::selectInput("enum_plot_group", "Group", choices = NULL),
            shiny::selectInput("enum_plot_k", "Profiles", choices = NULL),
            shiny::selectInput("enum_plot_type", "Type", c("bar", "line", "radar")),
            shiny::selectInput("enum_plot_palette", "Palette", names(mglpa_palettes))
          ),
          plot_box("enum_profile_plot", "480px"),
          shiny::downloadButton("dl_enum_plot_png", "Profile plot (PNG)", class = "btn-sm")
        ),
        table_card("Classification quality of the selected solution",
                   shiny::uiOutput("enum_classif_text"), DT::DTOutput("enum_classif"))
      )
    ),
    bslib::nav_panel(
      "3. Similarity tests",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::numericInput("sim_k", "Number of profiles retained", 3, min = 2, max = 10),
          shiny::selectInput("sim_var", "Within-profile variances",
                             c("Equal across profiles (class-invariant)" = "class_invariant",
                               "Profile-specific (class-varying)" = "class_varying")),
          starts_inputs("sim"),
          shiny::actionButton("run_sim", "Run similarity tests", class = "btn-primary w-100"),
          shiny::hr(),
          shiny::downloadButton("dl_sim_csv", "Table (CSV)", class = "btn-sm"),
          shiny::downloadButton("dl_sim_doc", "Table (Word)", class = "btn-sm")
        ),
        table_card("Fit indices and decisions", DT::DTOutput("sim_table"), shiny::uiOutput("sim_decision")),
        shiny::uiOutput("sim_notes"),
        bslib::card(
          bslib::card_header("Profile distribution and plot"),
          full_screen = TRUE,
          bslib::layout_columns(
            shiny::selectInput("sim_plot_model", "Model", choices = NULL),
            shiny::selectInput("sim_plot_type", "Type", c("bar", "line", "radar")),
            shiny::selectInput("sim_plot_palette", "Palette", names(mglpa_palettes))
          ),
          shiny::tableOutput("sim_props"),
          plot_box("sim_profile_plot", "480px"),
          shiny::downloadButton("dl_sim_plot_png", "Profile plot (PNG)", class = "btn-sm"),
          shiny::downloadButton("dl_sim_means_csv", "Profile means (CSV)", class = "btn-sm")
        ),
        table_card("Classification quality of the selected model",
                   shiny::uiOutput("sim_classif_text"), DT::DTOutput("sim_classif"))
      )
    ),
    bslib::nav_panel(
      "4. Predictive",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::selectInput("pred_source", "Profile solution held fixed", choices = NULL),
          shiny::uiOutput("pred_predictors_ui"),
          starts_inputs("pred", 20, 5),
          shiny::actionButton("run_pred", "Run predictive similarity", class = "btn-primary w-100"),
          shiny::hr(),
          shiny::downloadButton("dl_pred_csv", "Effects (CSV)", class = "btn-sm"),
          shiny::downloadButton("dl_pred_doc", "Effects (Word)", class = "btn-sm")
        ),
        table_card("Fit indices and decision", DT::DTOutput("pred_table"), shiny::uiOutput("pred_decision")),
        shiny::uiOutput("pred_notes"),
        table_card("Group-specific effects", DT::DTOutput("pred_effects_free")),
        table_card("Effects constrained equal across groups", DT::DTOutput("pred_effects_equal"))
      )
    ),
    bslib::nav_panel(
      "5. Explanatory",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::selectInput("expl_source", "Profile solution held fixed", choices = NULL),
          shiny::uiOutput("expl_outcomes_ui"),
          shiny::selectInput("expl_var", "Residual variances of the outcomes",
                             c("Both specifications" = "both", "One variance per outcome" = "equal",
                               "Profile- and group-specific" = "free")),
          starts_inputs("expl", 20, 5),
          shiny::actionButton("run_expl", "Run explanatory similarity", class = "btn-primary w-100"),
          shiny::hr(),
          shiny::downloadButton("dl_expl_csv", "Contrasts (CSV)", class = "btn-sm"),
          shiny::downloadButton("dl_expl_doc", "Contrasts (Word)", class = "btn-sm")
        ),
        table_card("Fit indices and decisions", DT::DTOutput("expl_table"), shiny::uiOutput("expl_decision")),
        shiny::uiOutput("expl_notes"),
        table_card("Outcome means and pairwise differences",
                   shiny::selectInput("expl_model", "Model", choices = NULL),
                   DT::DTOutput("expl_means"), DT::DTOutput("expl_contrasts"))
      )
    ),
    bslib::nav_panel(
      "6. Project",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 360,
          shiny::downloadButton("dl_project", "Save project (.rds)"),
          shiny::fileInput("load_project", "Load project (.rds)", accept = ".rds"),
          shiny::hr(),
          shiny::downloadButton("dl_report", "All tables (Word)")
        ),
        bslib::card(bslib::card_header("Project status"), shiny::uiOutput("project_status"))
      )
    )
  )
}

# ------------------------------------------------------------------------------
# Server -----------------------------------------------------------------------
# ------------------------------------------------------------------------------

app_server <- function(initial_data = NULL) {
  function(input, output, session) {
    rv <- shiny::reactiveValues(data = initial_data, enum = NULL, sim = NULL, pred = NULL,
                                expl = NULL, log = character())

    notify_error <- function(e) shiny::showNotification(paste("Error:", conditionMessage(e)),
                                                        type = "error", duration = 10)

    run_with_progress <- function(expr, message, n_steps) {
      if (inherits(session, "MockShinySession")) return(suppressMessages(expr))
      shiny::withProgress(message = message, value = 0, {
        step <- 0
        withCallingHandlers(expr, message = function(m) {
          step <<- step + 1
          shiny::setProgress(value = min(step / n_steps, 0.99), detail = trimws(conditionMessage(m)))
          invokeRestart("muffleMessage")
        })
      })
    }

    starts_of <- function(id) list(
      starts = c(input[[paste0(id, "_starts1")]], input[[paste0(id, "_starts2")]]),
      stiterations = input[[paste0(id, "_stiter")]], seed = input[[paste0(id, "_seed")]])
    n_cores <- function() { n <- input$cores; if (is.null(n) || is.na(n)) 1L else max(1L, as.integer(n)) }
    # Never hand a zero-width device to renderPlot (see plot_box in app_ui()).
    plot_width <- function(id) function() {
      w <- session$clientData[[paste0("output_", id, "_width")]]
      if (is.null(w) || is.na(w) || w < 100) 900 else w
    }

    # Fit-index table with convergence and replication badges. Replication of
    # the best log-likelihood is only meaningful for models estimated with
    # random starts (enumeration and similarity tests); the predictive and
    # explanatory models start from the retained solution, as with STARTS = 0
    # in Mplus, so only their convergence is reported (through the notes).
    fit_table_dt <- function(tab, badges = TRUE) {
      num <- intersect(c("LL", "AIC", "BIC", "CAIC", "ABIC", "entropy", "LRT", "VLMR_p", "aLMR_p",
                         "BLRT_p", "dCAIC", "dBIC", "dABIC"), names(tab))
      tab[num] <- lapply(tab[num], round, 3)
      if ("min_proportion" %in% names(tab)) tab$min_proportion <- round(tab$min_proportion, 3)
      tab$reps_n <- tab$replications
      tab$replicated <- paste0(tab$replications, "/", tab$n_final)
      tab$converged <- ifelse(tab$converged, "yes", "NO")
      cols <- intersect(c("group", "model", "k", "LL", "npar", "AIC", "BIC", "CAIC", "ABIC",
                          "entropy", "LRT", "VLMR_p", "aLMR_p", "BLRT_p", "min_proportion",
                          "dCAIC", "dBIC", "dABIC", "n_decrease",
                          "decision", if (badges) c("replicated", "converged", "reps_n")), names(tab))
      tab <- tab[, cols, drop = FALSE]
      hidden <- which(cols == "reps_n") - 1L
      dt <- DT::datatable(tab, rownames = FALSE,
                          options = list(dom = "t", scrollX = TRUE, pageLength = 50,
                                         columnDefs = if (length(hidden))
                                           list(list(visible = FALSE, targets = hidden)) else list()))
      if (badges) {
        dt <- DT::formatStyle(dt, "converged", fontWeight = "bold",
                              backgroundColor = DT::styleEqual(c("yes", "NO"), c("#d1fae5", "#fecaca")))
        dt <- DT::formatStyle(dt, "replicated", valueColumns = "reps_n", fontWeight = "bold",
                              backgroundColor = DT::styleInterval(c(1, 2), c("#fecaca", "#fef3c7", "#d1fae5")))
      }
      dt
    }
    # Mplus-style estimation notes for a named list of fits
    # The notes are shown in their own card below the table (never inside the
    # table card, where they could overlap the rows); identical messages are
    # merged into one entry that lists the models concerned.
    estimation_notes_ui <- function(fits, labels = names(fits), ok) {
      has <- vapply(fits, function(f) length(f$warnings) > 0L, logical(1))
      if (!any(has))
        return(bslib::card(
          class = "border-success",
          bslib::card_header(style = "background-color: #d1e7dd;", "Estimation notes"),
          bslib::card_body(fillable = FALSE, shiny::p(class = "mb-0", ok))))
      msgs <- unlist(lapply(which(has), function(i) fits[[i]]$warnings))
      labs <- unlist(lapply(which(has), function(i) rep(labels[i], length(fits[[i]]$warnings))))
      items <- lapply(unique(msgs), function(m)
        shiny::tags$li(class = "mb-1", shiny::strong(paste(labs[msgs == m], collapse = "; ")), ": ", m))
      n <- length(msgs)
      bslib::card(
        class = "border-warning",
        bslib::card_header(style = "background-color: #fff3cd;", "Estimation notes",
                           shiny::span(class = "badge text-bg-warning ms-2",
                                       sprintf("%d warning%s", n, if (n > 1L) "s" else ""))),
        bslib::card_body(fillable = FALSE, shiny::tags$ul(class = "mb-0 ps-3", items)))
    }

    # ---- Data ----------------------------------------------------------------
    shiny::observeEvent(input$file, {
      tryCatch({
        first <- readLines(input$file$datapath, n = 1, warn = FALSE)
        sep <- if (grepl(";", first, fixed = TRUE)) ";" else if (grepl("\t", first, fixed = TRUE)) "\t" else ","
        d <- utils::read.csv(input$file$datapath, sep = sep, dec = input$dec, check.names = FALSE,
                             strip.white = TRUE)
        if (nzchar(input$missing_code)) {
          code <- suppressWarnings(as.numeric(input$missing_code))
          for (v in names(d)) if (is.numeric(d[[v]]) && !is.na(code)) d[[v]][d[[v]] == code] <- NA
        }
        rv$data <- d
        rv$enum <- rv$sim <- rv$pred <- rv$expl <- NULL
      }, error = notify_error)
    })
    shiny::observeEvent(input$use_example, {
      e <- new.env()
      utils::data("mglpa_wm", package = "mglpa", envir = e)
      rv$data <- e$mglpa_wm
      rv$enum <- rv$sim <- rv$pred <- rv$expl <- NULL
      shiny::updateSelectizeInput(session, "indicators",
        selected = c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed"))
      shiny::updateSelectizeInput(session, "group", selected = "group")
      shiny::updateSelectizeInput(session, "predictors", selected = c("education", "cog_reserve"))
      shiny::updateSelectizeInput(session, "outcomes", selected = c("task_perf", "everyday_func"))
    }, ignoreInit = TRUE)
    shiny::observeEvent(rv$data, {
      d <- rv$data
      num <- names(d)[vapply(d, is.numeric, logical(1))]
      keep <- function(id) intersect(shiny::isolate(input[[id]]), names(d))
      shiny::updateSelectizeInput(session, "indicators", choices = num, selected = keep("indicators"))
      shiny::updateSelectizeInput(session, "group", choices = c("(none)" = "", names(d)),
                                  selected = keep("group"))
      shiny::updateSelectizeInput(session, "predictors", choices = num, selected = keep("predictors"))
      shiny::updateSelectizeInput(session, "outcomes", choices = num, selected = keep("outcomes"))
    })
    output$data_summary <- shiny::renderUI({
      d <- rv$data
      if (is.null(d)) return(shiny::p("No data loaded. Upload a CSV file or use the example data set."))
      g <- input$group
      txt <- sprintf("%d observations, %d variables.", nrow(d), ncol(d))
      if (!is.null(g) && nzchar(g) && g %in% names(d)) {
        tb <- table(d[[g]], useNA = "no")
        txt <- paste(txt, sprintf("Groups (%s): %s.", g, paste(sprintf("%s = %d", names(tb), tb), collapse = ", ")))
      }
      shiny::p(txt)
    })
    output$data_preview <- DT::renderDT({
      shiny::req(rv$data)
      DT::datatable(utils::head(rv$data, 200), options = list(pageLength = 8, scrollX = TRUE),
                    rownames = FALSE)
    })
    group_var <- shiny::reactive({ g <- input$group; if (is.null(g) || !nzchar(g)) NULL else g })

    # ---- Enumeration ---------------------------------------------------------
    shiny::observeEvent(input$run_enum, {
      shiny::req(rv$data, length(input$indicators) > 0)
      st <- starts_of("enum")
      ks <- seq(input$enum_k[1], input$enum_k[2])
      ng <- if (is.null(group_var())) 1 else length(unique(stats::na.omit(rv$data[[group_var()]])))
      tryCatch({
        rv$enum <- run_with_progress(
          enumerate_profiles(rv$data, input$indicators, group = group_var(), k = ks,
                             variances = input$enum_var, starts = st$starts,
                             stiterations = st$stiterations, seed = st$seed, verbose = TRUE,
                             cores = n_cores()),
          "Class enumeration", length(ks) * ng)
        shiny::updateSelectInput(session, "enum_plot_group", choices = names(rv$enum$fits))
        shiny::updateSelectInput(session, "enum_plot_k", choices = as.character(ks),
                                 selected = as.character(ks[min(3, length(ks))]))
        shiny::showNotification("Class enumeration finished.", type = "message")
      }, error = notify_error)
    })
    output$enum_table <- DT::renderDT({ shiny::req(rv$enum); fit_table_dt(rv$enum$table) })
    output$enum_notes <- shiny::renderUI({
      shiny::req(rv$enum)
      fits <- unlist(rv$enum$fits, recursive = FALSE)
      labels <- unlist(lapply(names(rv$enum$fits), function(g)
        paste0(if (g == "all") "" else paste0("Group ", g, ", "), "k = ", names(rv$enum$fits[[g]]))))
      # the BLRT objects carry their own local-maxima warnings
      if (length(rv$enum$blrt)) {
        fits <- c(fits, rv$enum$blrt)
        labels <- c(labels, vapply(names(rv$enum$blrt), function(nm) {
          parts <- strsplit(nm, " ", fixed = TRUE)[[1]]
          k <- as.integer(parts[length(parts)])
          g <- paste(parts[-length(parts)], collapse = " ")
          paste0("BLRT ", if (g == "all") "" else paste0("group ", g, ", "), k, " vs ", k - 1L, " profiles")
        }, character(1)))
      }
      estimation_notes_ui(fits, labels,
        ok = "All solutions converged and replicated their best log-likelihood at least three times.")
    })
    output$elbow <- shiny::renderPlot({ shiny::req(rv$enum); elbow_plot(rv$enum) },
                                      width = plot_width("elbow"), height = 420, res = 96)
    enum_fit <- shiny::reactive({
      shiny::req(rv$enum, input$enum_plot_group, input$enum_plot_k)
      rv$enum$fits[[input$enum_plot_group]][[input$enum_plot_k]]
    })
    output$enum_profile_plot <- shiny::renderPlot({
      f <- enum_fit(); shiny::req(f)
      profile_plot(f, type = input$enum_plot_type, palette = input$enum_plot_palette,
                   title = sprintf("%s: %s-profile solution",
                                   if (input$enum_plot_group == "all") "All" else paste("Group", input$enum_plot_group),
                                   input$enum_plot_k))
    }, width = plot_width("enum_profile_plot"), height = 480, res = 96)
    output$dl_enum_csv <- shiny::downloadHandler("enumeration.csv", function(file)
      utils::write.csv(rv$enum$table, file, row.names = FALSE))
    output$dl_enum_doc <- shiny::downloadHandler("enumeration.doc", function(file)
      write_table_doc(apa_table(rv$enum), file, title = "Table 1. Fit indices of the class enumeration"))
    output$dl_elbow_png <- shiny::downloadHandler("elbow_plot.png", function(file)
      ggplot2::ggsave(file, elbow_plot(rv$enum), width = 9, height = 4.5, dpi = 300))
    output$dl_enum_plot_png <- shiny::downloadHandler("profile_plot.png", function(file)
      ggplot2::ggsave(file, profile_plot(enum_fit(), type = input$enum_plot_type,
                                         palette = input$enum_plot_palette), width = 9, height = 5, dpi = 300))
    # Bootstrap LRT (TECH14) for all consecutive solutions of the enumeration
    shiny::observeEvent(input$run_blrt, {
      shiny::req(rv$enum)
      n_comp <- sum(vapply(rv$enum$fits, function(fg) sum(diff(as.integer(names(fg))) == 1L),
                           integer(1)))
      if (n_comp == 0L) {
        shiny::showNotification("The BLRT needs at least two consecutive solutions.", type = "warning")
        return()
      }
      B <- max(1L, as.integer(input$blrt_B))
      st <- c(max(1L, as.integer(input$blrt_starts1)), max(1L, as.integer(input$blrt_starts2)))
      seed <- input$enum_seed
      if (is.null(seed) || is.na(seed)) seed <- NULL
      tryCatch({
        rv$enum <- run_with_progress(
          enumeration_blrt(rv$enum, B = B, starts = st, seed = seed, verbose = TRUE,
                           cores = n_cores()),
          "Bootstrap LRT", n_comp * ceiling(B / blrt_chunk(n_cores())))
        shiny::showNotification("Bootstrap LRT finished.", type = "message")
      }, error = notify_error)
    })
    # Classification quality: entropy and average posterior probabilities
    classif_dt <- function(fit) {
      tab <- apa_table(classification_quality(fit))
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 50))
    }
    classif_text <- function(fit) {
      cq <- classification_quality(fit)
      if (is.na(cq$entropy)) return(shiny::p("Entropy is not defined for a single profile."))
      shiny::p(shiny::strong(sprintf("Entropy = %.3f", cq$entropy)),
               if (length(cq$group_levels) > 1L)
                 sprintf(" (profiles x groups; profile-specific entropy = %.3f; by group: %s)",
                         cq$entropy_profiles,
                         paste(sprintf("%s = %.3f", names(cq$entropy_by_group), cq$entropy_by_group),
                               collapse = ", "))
               else "",
               ". Rows: observations assigned to each profile; columns: their average posterior",
               " probabilities of membership in each profile.")
    }
    output$enum_classif_text <- shiny::renderUI({ f <- enum_fit(); shiny::req(f); classif_text(f) })
    output$enum_classif <- DT::renderDT({ f <- enum_fit(); shiny::req(f); classif_dt(f) })

    # ---- Similarity ----------------------------------------------------------
    shiny::observeEvent(input$run_sim, {
      shiny::req(rv$data, length(input$indicators) > 0)
      if (is.null(group_var())) {
        shiny::showNotification("Select a grouping variable in the Data tab.", type = "warning"); return()
      }
      st <- starts_of("sim")
      tryCatch({
        rv$sim <- run_with_progress(
          similarity_tests(rv$data, input$indicators, group = group_var(), k = input$sim_k,
                           variances = input$sim_var, starts = st$starts,
                           stiterations = st$stiterations, seed = st$seed, verbose = TRUE,
                           cores = n_cores()),
          "Similarity tests", 4)
        rv$pred <- rv$expl <- NULL
        nms <- names(rv$sim$fits)
        shiny::updateSelectInput(session, "sim_plot_model", choices = nms, selected = rv$sim$retained)
        shiny::updateSelectInput(session, "pred_source", choices = nms, selected = rv$sim$retained)
        shiny::updateSelectInput(session, "expl_source", choices = nms, selected = rv$sim$retained)
        shiny::showNotification("Similarity tests finished.", type = "message")
      }, error = notify_error)
    })
    output$sim_table <- DT::renderDT({ shiny::req(rv$sim); fit_table_dt(rv$sim$table) })
    output$sim_notes <- shiny::renderUI({
      shiny::req(rv$sim)
      estimation_notes_ui(rv$sim$fits,
        ok = "All models converged and replicated their best log-likelihood at least three times.")
    })
    output$sim_decision <- shiny::renderUI({
      shiny::req(rv$sim)
      shiny::p(shiny::strong("Most constrained model supported: "), rv$sim$retained,
               ". Decision rule: at least two of CAIC, BIC and ABIC lower than the previous model.")
    })
    sim_fit <- shiny::reactive({ shiny::req(rv$sim, input$sim_plot_model); rv$sim$fits[[input$sim_plot_model]] })
    output$sim_props <- shiny::renderTable({
      f <- sim_fit(); shiny::req(f)
      pr <- profile_proportions(f)
      out <- data.frame(Group = rownames(pr), round(100 * pr, 1), check.names = FALSE)
      names(out)[-1] <- paste0("Profile ", seq_len(ncol(pr)), " (%)")
      out
    })
    output$sim_classif_text <- shiny::renderUI({ f <- sim_fit(); shiny::req(f); classif_text(f) })
    output$sim_classif <- DT::renderDT({ f <- sim_fit(); shiny::req(f); classif_dt(f) })
    output$sim_profile_plot <- shiny::renderPlot({
      f <- sim_fit(); shiny::req(f)
      profile_plot(f, type = input$sim_plot_type, palette = input$sim_plot_palette,
                   title = paste(input$sim_plot_model, "similarity model"))
    }, width = plot_width("sim_profile_plot"), height = 480, res = 96)
    output$dl_sim_csv <- shiny::downloadHandler("similarity_tests.csv", function(file)
      utils::write.csv(rv$sim$table, file, row.names = FALSE))
    output$dl_sim_doc <- shiny::downloadHandler("similarity_tests.doc", function(file)
      write_table_doc(apa_table(rv$sim), file, title = "Table 2. Fit indices of the multi-group similarity tests"))
    output$dl_sim_plot_png <- shiny::downloadHandler("profile_plot.png", function(file)
      ggplot2::ggsave(file, profile_plot(sim_fit(), type = input$sim_plot_type,
                                         palette = input$sim_plot_palette), width = 9, height = 5, dpi = 300))
    output$dl_sim_means_csv <- shiny::downloadHandler("profile_means.csv", function(file)
      utils::write.csv(profile_means(sim_fit()), file, row.names = FALSE))

    # ---- Predictive ----------------------------------------------------------
    output$pred_predictors_ui <- shiny::renderUI({
      shiny::req(rv$data)
      num <- names(rv$data)[vapply(rv$data, is.numeric, logical(1))]
      shiny::selectizeInput("pred_vars", "Predictors", choices = num, selected = input$predictors,
                            multiple = TRUE)
    })
    shiny::observeEvent(input$run_pred, {
      shiny::req(rv$sim, length(input$pred_vars) > 0)
      st <- starts_of("pred")
      tryCatch({
        rv$pred <- run_with_progress(
          predictive_similarity(rv$sim, predictors = input$pred_vars, source = input$pred_source,
                                starts = st$starts, stiterations = st$stiterations, seed = st$seed,
                                verbose = TRUE),
          "Predictive similarity", 2)
        shiny::showNotification("Predictive similarity finished.", type = "message")
      }, error = notify_error)
    })
    effects_dt <- function(e) {
      e$coef <- round(e$coef, 3); e$se <- round(e$se, 3); e$z <- round(e$z, 2); e$p <- round(e$p, 3)
      e$OR <- round(e$OR, 3); e$OR_ci_lower <- round(e$OR_ci_lower, 3); e$OR_ci_upper <- round(e$OR_ci_upper, 3)
      DT::datatable(e[, c("group_level", "predictor", "comparison", "coef", "se", "z", "p", "OR",
                          "OR_ci_lower", "OR_ci_upper")],
                    options = list(pageLength = 25, dom = "t", scrollX = TRUE), rownames = FALSE)
    }
    output$pred_table <- DT::renderDT({ shiny::req(rv$pred); fit_table_dt(rv$pred$table, badges = FALSE) })
    output$pred_notes <- shiny::renderUI({
      shiny::req(rv$pred)
      fits <- rv$pred$fits
      names(fits) <- paste(names(fits), "effects")
      estimation_notes_ui(fits, ok = "Both models converged from the retained profile solution.")
    })
    output$pred_decision <- shiny::renderUI({
      shiny::req(rv$pred, rv$pred$decision)
      shiny::p(shiny::strong("Predictive similarity: "), rv$pred$decision$decision,
               sprintf(" (%d of 3 information criteria lower for equal effects).", rv$pred$decision$n_decrease))
    })
    output$pred_effects_free <- DT::renderDT({ shiny::req(rv$pred$effects$free); effects_dt(rv$pred$effects$free) })
    output$pred_effects_equal <- DT::renderDT({ shiny::req(rv$pred$effects$equal); effects_dt(rv$pred$effects$equal) })
    output$dl_pred_csv <- shiny::downloadHandler("predictive_effects.csv", function(file) {
      e <- do.call(rbind, lapply(names(rv$pred$effects), function(nm) cbind(model = nm, rv$pred$effects[[nm]])))
      utils::write.csv(e, file, row.names = FALSE)
    })
    output$dl_pred_doc <- shiny::downloadHandler("predictive_effects.doc", function(file)
      write_table_doc(apa_table(rv$pred), file, title = "Table 3. Predictors of profile membership"))

    # ---- Explanatory ---------------------------------------------------------
    output$expl_outcomes_ui <- shiny::renderUI({
      shiny::req(rv$data)
      num <- names(rv$data)[vapply(rv$data, is.numeric, logical(1))]
      shiny::selectizeInput("expl_vars", "Outcomes", choices = num, selected = input$outcomes, multiple = TRUE)
    })
    shiny::observeEvent(input$run_expl, {
      shiny::req(rv$sim, length(input$expl_vars) > 0)
      st <- starts_of("expl")
      tryCatch({
        rv$expl <- run_with_progress(
          explanatory_similarity(rv$sim, outcomes = input$expl_vars, source = input$expl_source,
                                 outcome_variances = input$expl_var, starts = st$starts,
                                 stiterations = st$stiterations, seed = st$seed, verbose = TRUE),
          "Explanatory similarity", if (input$expl_var == "both") 4 else 2)
        shiny::updateSelectInput(session, "expl_model", choices = names(rv$expl$fits))
        shiny::showNotification("Explanatory similarity finished.", type = "message")
      }, error = notify_error)
    })
    output$expl_table <- DT::renderDT({ shiny::req(rv$expl); fit_table_dt(rv$expl$table, badges = FALSE) })
    output$expl_notes <- shiny::renderUI({
      shiny::req(rv$expl)
      estimation_notes_ui(rv$expl$fits, ok = "All models converged from the retained profile solution.")
    })
    output$expl_decision <- shiny::renderUI({
      shiny::req(rv$expl)
      shiny::tagList(lapply(names(rv$expl$decisions), function(ov)
        shiny::p(shiny::strong(sprintf("Explanatory similarity (%s outcome variances): ", ov)),
                 rv$expl$decisions[[ov]]$decision,
                 sprintf(" (%d of 3 criteria lower for equal means).", rv$expl$decisions[[ov]]$n_decrease))))
    })
    output$expl_means <- DT::renderDT({
      shiny::req(rv$expl, input$expl_model)
      m <- rv$expl$means[[input$expl_model]]
      m[c("mean", "se", "ci_lower", "ci_upper")] <- lapply(m[c("mean", "se", "ci_lower", "ci_upper")], round, 3)
      m$proportion <- round(m$proportion, 3)
      DT::datatable(m[, c("group_level", "class", "variable", "mean", "se", "ci_lower", "ci_upper", "proportion")],
                    options = list(dom = "t", pageLength = 50), rownames = FALSE, caption = "Outcome means by profile and group")
    })
    output$expl_contrasts <- DT::renderDT({
      shiny::req(rv$expl, input$expl_model)
      ct <- rv$expl$contrasts[[input$expl_model]]
      ct[c("estimate", "se", "ci_lower", "ci_upper")] <- lapply(ct[c("estimate", "se", "ci_lower", "ci_upper")], round, 3)
      ct$z <- round(ct$z, 2); ct$p <- round(ct$p, 3)
      DT::datatable(ct[, c("outcome", "type", "group_level", "comparison", "estimate", "se", "z", "p", "ci_lower", "ci_upper")],
                    options = list(dom = "t", pageLength = 100), rownames = FALSE, caption = "Pairwise differences (delta method)")
    })
    output$dl_expl_csv <- shiny::downloadHandler("explanatory_contrasts.csv", function(file) {
      e <- do.call(rbind, lapply(names(rv$expl$contrasts), function(nm) cbind(model = nm, rv$expl$contrasts[[nm]])))
      utils::write.csv(e, file, row.names = FALSE)
    })
    output$dl_expl_doc <- shiny::downloadHandler("explanatory_contrasts.doc", function(file)
      write_table_doc(apa_table(rv$expl), file, title = "Table 4. Differences in outcome means"))

    # ---- Project -------------------------------------------------------------
    output$dl_project <- shiny::downloadHandler(
      function() paste0("mglpa_project_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"),
      function(file) saveRDS(list(version = as.character(utils::packageVersion("mglpa")),
                                  saved = Sys.time(), data = rv$data,
                                  settings = list(indicators = input$indicators, group = input$group,
                                                  predictors = input$predictors, outcomes = input$outcomes),
                                  enum = rv$enum, sim = rv$sim, pred = rv$pred, expl = rv$expl), file))
    shiny::observeEvent(input$load_project, {
      tryCatch({
        pr <- readRDS(input$load_project$datapath)
        rv$data <- pr$data
        rv$enum <- pr$enum; rv$sim <- pr$sim; rv$pred <- pr$pred; rv$expl <- pr$expl
        s <- pr$settings
        shiny::updateSelectizeInput(session, "indicators", selected = s$indicators)
        shiny::updateSelectizeInput(session, "group", selected = s$group)
        shiny::updateSelectizeInput(session, "predictors", selected = s$predictors)
        shiny::updateSelectizeInput(session, "outcomes", selected = s$outcomes)
        if (!is.null(pr$enum)) {
          shiny::updateSelectInput(session, "enum_plot_group", choices = names(pr$enum$fits))
          shiny::updateSelectInput(session, "enum_plot_k", choices = as.character(pr$enum$k))
        }
        if (!is.null(pr$sim)) {
          nms <- names(pr$sim$fits)
          shiny::updateSelectInput(session, "sim_plot_model", choices = nms, selected = pr$sim$retained)
          shiny::updateSelectInput(session, "pred_source", choices = nms, selected = pr$sim$retained)
          shiny::updateSelectInput(session, "expl_source", choices = nms, selected = pr$sim$retained)
        }
        if (!is.null(pr$expl)) shiny::updateSelectInput(session, "expl_model", choices = names(pr$expl$fits))
        shiny::showNotification("Project loaded.", type = "message")
      }, error = notify_error)
    })
    output$project_status <- shiny::renderUI({
      items <- c(
        sprintf("Data: %s", if (is.null(rv$data)) "none" else sprintf("%d x %d", nrow(rv$data), ncol(rv$data))),
        sprintf("Class enumeration: %s", if (is.null(rv$enum)) "not run" else sprintf("%d models", nrow(rv$enum$table))),
        sprintf("Similarity tests: %s", if (is.null(rv$sim)) "not run" else sprintf("k = %d, retained: %s", rv$sim$k, rv$sim$retained)),
        sprintf("Predictive similarity: %s", if (is.null(rv$pred)) "not run" else rv$pred$decision$decision),
        sprintf("Explanatory similarity: %s", if (is.null(rv$expl)) "not run" else
          paste(vapply(names(rv$expl$decisions), function(ov) sprintf("%s variances: %s", ov, rv$expl$decisions[[ov]]$decision), ""), collapse = "; ")))
      shiny::tags$ul(lapply(items, shiny::tags$li))
    })
    output$dl_report <- shiny::downloadHandler("mglpa_tables.doc", function(file) {
      parts <- list()
      if (!is.null(rv$enum)) parts[["Table 1. Fit indices of the class enumeration"]] <- apa_table(rv$enum)
      if (!is.null(rv$sim)) parts[["Table 2. Fit indices of the multi-group similarity tests"]] <- apa_table(rv$sim)
      if (!is.null(rv$sim))
        parts[[sprintf("Table 2b. Classification quality of the %s similarity model", rv$sim$retained)]] <-
          apa_table(classification_quality(rv$sim$fits[[rv$sim$retained]]))
      if (!is.null(rv$pred)) parts[["Table 3. Predictors of profile membership"]] <- apa_table(rv$pred)
      if (!is.null(rv$expl)) parts[["Table 4. Differences in outcome means"]] <- apa_table(rv$expl)
      if (!length(parts)) stop("No results to export.")
      write_tables_doc(parts, file)
    })
  }
}
