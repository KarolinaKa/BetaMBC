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
  multi_beta_clus(NonBasalLike, groups = 3, location, scale )


table("BasalLike" = cluster_BasalLike$ClusterMemberships, "NonBasalLike" = cluster_NonBasalLike$ClusterMemberships)
1 - sum(diag(tab))/sum(tab)

plot(cluster_BasalLike, freq = T)
plot(cluster_NonBasalLike, freq = T)


ranks <- ranking(cluster_BasalLike, cluster_NonBasalLike, RowNames = row.names(BasalLike))
plot(ranks, Rank = 1, BasalLike, NonBasalLike)

