ind <- c("verbal_wm", "visual_wm", "updating", "inhibition", "proc_speed")
wm <- mglpa::mglpa_wm
wm1 <- wm[wm$group == 1, ]

# Reference values parsed from the Mplus outputs of the example project
mplus_ref <- read.csv(system.file("extdata", "mplus_reference_fits.csv", package = "mglpa"),
                      stringsAsFactors = FALSE)
# Mplus TECH11 (VLMR/aLMR) and TECH14 (BLRT) values of the enumeration runs
mplus_t11 <- read.csv(system.file("extdata", "mplus_tech11.csv", package = "mglpa"),
                      stringsAsFactors = FALSE, colClasses = c(mplus_version = "character"))
ref_value <- function(step, column, group = NA, k = NA, model = NA) {
  r <- mplus_ref[mplus_ref$step == step, ]
  if (!is.na(group)) r <- r[r$group == group, ]
  if (!is.na(k)) r <- r[r$k == k, ]
  if (!is.na(model)) r <- r[r$model == model, ]
  r[[column]][1]
}
