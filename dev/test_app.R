# Headless check of the Shiny application with shiny::testServer.
pkgload::load_all("mglpa", quiet = TRUE)
library(shiny)
ui <- mglpa:::app_ui()
stopifnot(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
cat("UI built\n")
server <- mglpa:::app_server(NULL)
testServer(server, {
  session$setInputs(dec = ".", missing_code = "")
  session$setInputs(use_example = 1)
  stopifnot(!is.null(rv$data))
  session$setInputs(indicators = c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed"),
                    group = "group", predictors = c("education", "cog_reserve"),
                    outcomes = c("task_perf", "everyday_func"))
  cat("data ok:", nrow(rv$data), "rows\n")
  stopifnot(grepl("Groups", as.character(output$data_summary$html)))

  session$setInputs(enum_k = c(1, 3), enum_var = "class_invariant", enum_starts1 = 10,
                    enum_starts2 = 2, enum_stiter = 10, enum_seed = 1)
  session$setInputs(run_enum = 1)
  stopifnot(!is.null(rv$enum), nrow(rv$enum$table) == 6)
  cat("enumeration ok\n")
  session$setInputs(enum_plot_group = "1", enum_plot_k = "3", enum_plot_type = "bar",
                    enum_plot_palette = "vibrant")
  stopifnot(!is.null(output$enum_profile_plot))
  stopifnot(!is.null(output$elbow))
  stopifnot(!is.null(output$enum_table))
  stopifnot(all(c("VLMR_p", "aLMR_p") %in% names(rv$enum$table)))
  stopifnot(!is.null(output$enum_classif), !is.null(output$enum_classif_text))
  cat("enumeration outputs ok\n")
  session$setInputs(blrt_B = 3, blrt_starts1 = 3, blrt_starts2 = 1)
  session$setInputs(run_blrt = 1)
  stopifnot("BLRT_p" %in% names(rv$enum$table), !anyNA(rv$enum$table$BLRT_p[rv$enum$table$k > 1]))
  cat("blrt ok\n")

  session$setInputs(sim_k = 3, sim_var = "class_invariant", sim_starts1 = 20, sim_starts2 = 4,
                    sim_stiter = 10, sim_seed = 2)
  session$setInputs(run_sim = 1)
  stopifnot(!is.null(rv$sim), rv$sim$retained == "dispersion")
  cat("similarity ok:", rv$sim$retained, "\n")
  session$setInputs(sim_plot_model = "dispersion", sim_plot_type = "line", sim_plot_palette = "apa")
  stopifnot(!is.null(output$sim_profile_plot), !is.null(output$sim_props), !is.null(output$sim_table))
  stopifnot(!is.null(output$sim_classif), !is.null(output$sim_classif_text))
  stopifnot(grepl("dispersion", as.character(output$sim_decision$html)))

  session$setInputs(pred_source = "dispersion", pred_vars = c("education", "cog_reserve"),
                    pred_starts1 = 5, pred_starts2 = 2, pred_stiter = 10, pred_seed = 3)
  session$setInputs(run_pred = 1)
  stopifnot(!is.null(rv$pred))
  cat("predictive ok:", rv$pred$decision$decision, "\n")
  stopifnot(!is.null(output$pred_table), !is.null(output$pred_effects_free), !is.null(output$pred_effects_equal))

  session$setInputs(expl_source = "dispersion", expl_vars = c("task_perf", "everyday_func"),
                    expl_var = "equal", expl_starts1 = 5, expl_starts2 = 2, expl_stiter = 10, expl_seed = 4)
  session$setInputs(run_expl = 1)
  stopifnot(!is.null(rv$expl))
  cat("explanatory ok:", rv$expl$decisions$equal$decision, "\n")
  session$setInputs(expl_model = "free_means_equal_var")
  stopifnot(!is.null(output$expl_means), !is.null(output$expl_contrasts), !is.null(output$expl_table))
  stopifnot(grepl("Explanatory", as.character(output$expl_decision$html)))
  stopifnot(grepl("retained", as.character(output$project_status$html)))

  tmp <- tempfile(fileext = ".doc")
  mglpa:::write_tables_doc(list("T1" = apa_table(rv$enum), "T2" = apa_table(rv$sim),
                                "T3" = apa_table(rv$pred), "T4" = apa_table(rv$expl)), tmp)
  stopifnot(file.exists(tmp), file.size(tmp) > 1000)
  cat("report ok:", file.size(tmp), "bytes\n")
  tmp2 <- tempfile(fileext = ".rds")
  saveRDS(list(data = rv$data, enum = rv$enum, sim = rv$sim, pred = rv$pred, expl = rv$expl,
               settings = list(indicators = input$indicators, group = input$group)), tmp2)
  session$setInputs(load_project = list(datapath = tmp2, name = "p.rds"))
  stopifnot(!is.null(rv$sim))
  cat("project load ok\n")
})
cat("APP TEST PASSED\n")
