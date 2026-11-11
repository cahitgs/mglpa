# End-to-end check of the high-level workflow on the example data.
pkgload::load_all("mglpa", quiet = TRUE)
data(mglpa_wm)
ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
dir.create("mglpa/dev/out", showWarnings = FALSE)

cat("\n==== Step 0: enumeration (k = 1:4, small starts) ====\n")
enum <- enumerate_profiles(mglpa_wm, ind, group = "group", k = 1:4, starts = c(30, 5), seed = 1)
print(enum)
p <- elbow_plot(enum)
ggplot2::ggsave("mglpa/dev/out/elbow.png", p, width = 8, height = 4, dpi = 100)
f3 <- enum$fits[["1"]][["3"]]
print(profile_means(f3))
ggplot2::ggsave("mglpa/dev/out/profile_bar.png", profile_plot(f3), width = 7, height = 4, dpi = 100)
ggplot2::ggsave("mglpa/dev/out/profile_line.png", profile_plot(f3, "line"), width = 7, height = 4, dpi = 100)
ggplot2::ggsave("mglpa/dev/out/profile_radar.png", profile_plot(f3, "radar"), width = 6, height = 5, dpi = 100)
print(apa_table(enum))

cat("\n==== Step 1: similarity tests (k = 3) ====\n")
sim <- similarity_tests(mglpa_wm, ind, group = "group", k = 3, starts = c(60, 10), seed = 2)
print(sim)
print(apa_table(sim))
write_table_doc(apa_table(sim), "mglpa/dev/out/similarity.doc", title = "Table 2. Fit indices of the similarity tests")
ggplot2::ggsave("mglpa/dev/out/profile_sim.png", profile_plot(sim$fits$dispersion), width = 9, height = 4, dpi = 100)
print(profile_proportions(sim$fits$dispersion))

cat("\n==== Step 2: predictive similarity ====\n")
pred <- predictive_similarity(sim, predictors = c("education", "cog_reserve"), seed = 3)
print(pred)
print(apa_table(pred))

cat("\n==== Step 3: explanatory similarity ====\n")
expl <- explanatory_similarity(sim, outcomes = c("task_perf", "everyday_func"), seed = 4)
print(expl)
print(expl$contrasts$equal_means_equal_var)
print(apa_table(expl))
write_table_doc(apa_table(expl), "mglpa/dev/out/explanatory.doc", title = "Table 4. Explanatory similarity")

cat("\n==== reorder check ====\n")
f <- sim$fits$dispersion
f2 <- reorder_profiles(f, c(3, 1, 2))
cat("LL before/after:", f$loglik, f2$loglik, "\n")
print(round(profile_proportions(f), 3)); print(round(profile_proportions(f2), 3))
cat("\nDone.\n")
