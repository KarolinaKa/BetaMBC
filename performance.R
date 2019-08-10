#################################################
## multi_beta_clus performance                 ##
##                                             ##
## date: July 2019                             ##
## Dataset file: "analysis_data.rds"           ##
## Supplementary file: "analysis_subtypes.rds" ##
#################################################

source("multi_beta_clus.R")
library(microbenchmark)

analysis_data <- readRDS("analysis_data.rds")
analysis_subtypes <- readRDS("analysis_subtypes.rds")

basal_like <-
  analysis_data[, analysis_subtypes == "Basal-like"]
non_basal_like <-
  analysis_data[, analysis_subtypes != "Basal-like"]

microbenchmark(multi_beta_clus(
  basal_like,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
), times = 20)

microbenchmark(multi_beta_clus(
  non_basal_like,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
), times = 20)

microbenchmark(multi_beta_clus(
  analysis_data,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
), times = 20)


