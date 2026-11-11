# Parallel random starts need the installed package (workers load the namespace).
# Run after: devtools::install("mglpa")
library(mglpa)
data(mglpa_wm)
ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
d1 <- mglpa_wm[mglpa_wm$group == 1, ]
cat("cores available:", parallel::detectCores(), "\n")
for (cores in c(1, 4)) {
  t0 <- Sys.time()
  f <- fit_lpa(d1, mglpa_model(8, ind), starts = c(3000, 300), seed = 18, cores = cores)
  cat(sprintf("k=8 3000/300 cores=%d: LL = %.3f  reps = %d/300  %.1fs\n", cores, f$loglik,
              f$starts$best_ll_replications, as.numeric(Sys.time() - t0, units = "secs")))
}
t0 <- Sys.time()
enum <- enumerate_profiles(mglpa_wm, ind, group = "group", k = 1:8, starts = c(1000, 100),
                           seed = 1, verbose = FALSE, cores = 4)
cat(sprintf("enumeration k=1:8 x 2 groups, 1000/100 starts, 4 cores: %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))
print(enum$table[, c("group", "k", "LL", "BIC", "replications")])
