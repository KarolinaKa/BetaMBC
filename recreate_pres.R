##############################################
## Recreate presentation results            ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: July 2019                          ##
##############################################

source("multi_beta_clus.R")
source("ranking.R")

BasalLike <- readRDS("FullBasalLike.rds")
NonBasalLike <- readRDS("FullNonBasalLike.rds")

location <- c(0, 0.5, 1)
scale <- c(0.1, 0.1, 0.1)

set.seed(13323)

cluster_BasalLike <-
  multi_beta_clus(BasalLike, groups = 3, location, scale)

cluster_NonBasalLike <-
  multi_beta_clus(NonBasalLike, groups = 3, location, scale)


table(
  "BasalLike" = cluster_BasalLike$ClusterMemberships,
  "NonBasalLike" = cluster_NonBasalLike$ClusterMemberships
)
1 - sum(diag(tab)) / sum(tab)

par(mfrow = c(1, 2))
plot(cluster_BasalLike, freq = T)
plot(cluster_NonBasalLike, freq = T)
par(mfrow = c(1, 1))

mean(cluster_BasalLike$ClusterUncertainties)
median(cluster_BasalLike$ClusterUncertainties)

mean(cluster_NonBasalLike$ClusterUncertainties)
median(cluster_NonBasalLike$ClusterUncertainties)

par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
hist(
  cluster_BasalLike$ClusterUncertainties,
  xlim = c(0, 0.7),
  xlab = "Uncertainty",
  main = "Cohort A"
)
hist(
  cluster_NonBasalLike$ClusterUncertainties,
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
  "Max. uncertainty is 0.6666",
  side = 1,
  outer = TRUE,
  cex = 1,
  font = 3
)

ranks <-
  ranking(cluster_BasalLike,
          cluster_NonBasalLike,
          RowNames = row.names(BasalLike))
plot(ranks, Rank = 1, BasalLike, NonBasalLike)

par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
hist(ranks[, 1],
     xlim = c(0, 0.7),
     xlab = "Uncertainty",
     main = "Cohort A")
hist(ranks[, 2],
     xlim = c(0, 0.7),
     xlab = "Uncertainty",
     main = "Cohort B")
mtext(
  "Uncertainty Histograms",
  outer = TRUE,
  cex = 1.3,
  font = 2
)
mtext(
  "Max. uncertainty is 0.6666",
  side = 1,
  outer = TRUE,
  cex = 1,
  font = 3
)
summary(ranks[, 1])
summary(ranks[, 2])



