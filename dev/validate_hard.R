pkgload::load_all("mglpa", quiet = TRUE)
d <- read.csv("mglpa/inst/extdata/mglpa_wm_dataset_n1000.csv")
ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
ref <- c("5" = -2161.878, "8" = -2141.698)
for (k in c(5, 8)) {
  m <- mglpa_model(k, ind)
  f <- fit_lpa(d[d$group == 1, ], m, starts = c(3000, 300), stiterations = 10, seed = 10 + k, se = "none")
  cat(sprintf("Group 1 k=%d  starts 3000/300: LL = %.3f (Mplus %.3f, diff %+.4f) reps = %d/%d  %.0fs\n",
              k, f$loglik, ref[as.character(k)], f$loglik - ref[as.character(k)],
              f$starts$best_ll_replications, length(f$starts$final_ll), f$time))
  cat("   top final LLs:", paste(round(head(f$starts$final_ll, 8), 3), collapse = " "), "\n")
  flush.console()
}
