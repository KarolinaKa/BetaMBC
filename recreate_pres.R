##############################################
## Recreate presentation results            ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: July 2019                          ##
##############################################

source("multi_beta_clus.R")
source("ranking.R")

cohortA <- readRDS("FullBasalLike.rds")
cohortB <- readRDS("FullNonBasalLike.rds")

location <- c(0, 0.5, 1)
scale <- c(0.1, 0.1, 0.1)

cluster_cohortA <-
  multi_beta_clus(cohortA, groups = 3, location, scale)

cluster_cohortB <-
  multi_beta_clus(cohortB, groups = 3, location, scale)

tab <- table(
  "cohortA" = cluster_cohortA$ClusterMemberships,
  "cohortB" = cluster_cohortB$ClusterMemberships
)
tab
1 - sum(diag(tab)) / sum(tab)

par(mfrow = c(1, 2))
plot(cluster_cohortA, freq = T)
plot(cluster_cohortB, freq = T)
par(mfrow = c(1, 1))

par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
hist(
  cluster_cohortA$ClusterUncertainties,
  xlim = c(0, 0.7),
  xlab = "Uncertainty",
  main = "Cohort A", 
  ylim = c(0, 12000)
)
hist(
  cluster_cohortB$ClusterUncertainties,
  xlim = c(0, 0.7),
  xlab = "Uncertainty",
  main = "Cohort B"
)
mtext(
  "Uncertainty Histograms",
  outer = TRUE,
  cex = 1.3,
  font = 2
)
mtext(
  "Max. possible uncertainty is 1 - 1/G = 0.67",
  side = 1,
  outer = TRUE,
  cex = 1,
  font = 3
)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

ranks <-
  ranking(cluster_cohortA,
          cluster_cohortB,
          RowNames = row.names(cohortA))
plot(ranks, Rank = 1, cohortA, cohortB)

par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
hist(ranks[, 1],
     xlim = c(0, 0.7),
     xlab = "Uncertainty",
     main = "Cohort A",
     ylim = c(0, 500))
hist(ranks[, 2],
     xlim = c(0, 0.7),
     xlab = "Uncertainty",
     main = "Cohort B",
     ylim = c(0, 500))
mtext(
  "Uncertainty Histograms",
  outer = TRUE,
  cex = 1.3,
  font = 2
)
mtext(
  "Max. possible uncertainty is 1 - 1/G = 0.67",
  side = 1,
  outer = TRUE,
  cex = 1,
  font = 3
)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

summary(ranks[, 1])
summary(ranks[, 2])



