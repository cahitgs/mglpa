#' Simulated working-memory data set (two age groups)
#'
#' The simulated data set distributed with the MG-LPA Analyzer: 1,000
#' participants (500 younger and 500 older adults) with five standardized
#' working-memory indicators, two predictors of profile membership and two
#' distal outcomes. The Mplus results reported for this data set in the
#' accompanying article are reproduced by the estimation functions of this
#' package.
#'
#' @format A data frame with 1,000 rows and 11 variables:
#' \describe{
#'   \item{id}{Participant identifier.}
#'   \item{group}{Age group: 1 = younger adults, 2 = older adults.}
#'   \item{verbal_wm, visual_wm, updating, inhibition, proc_speed}{Profile
#'     indicators (z scores).}
#'   \item{education}{Years of education (predictor).}
#'   \item{cog_reserve}{Cognitive reserve (z score; predictor).}
#'   \item{task_perf}{Task performance, percent correct (outcome).}
#'   \item{everyday_func}{Everyday functioning, 0-100 scale (outcome).}
#' }
#' @source Marangoz, M. C., Cokluk-Bokeoglu, O., & Morin, A. J. S. (2026).
#'   MG-LPA Analyzer (Version 1.0.0). \doi{10.5281/zenodo.22698717}
"mglpa_wm"
