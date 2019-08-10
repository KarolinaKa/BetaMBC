##############################################
## Setting optimal scale param              ##
## Multvariate Clustering                   ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

source('run_em.R')
source('multi_beta_clus.R')

analysis_data <- readRDS("analysis_data.rds")
analysis_subtypes <- readRDS("analysis_subtypes.rds")

basal_like <- analysis_data[, analysis_subtypes == "Basal-like"]
non_basal_like <- analysis_data[, analysis_subtypes != "Basal-like"]

# run MBC assuming small scale

location <- c(0, 0.5, 1)
scale <- c(0.01, 0.01, 0.01)

res_bl <- multi_beta_clus(basal_like, groups = 3, location, scale)
res_nbl <- multi_beta_clus(non_basal_like, groups = 3, location, scale)

# Basal-like group 

cluster1 <- which(res_bl$ClusterProbabilities[, 1] > 0.9)
cluster1 <- sample(cluster1, 500)
cluster1 <- basal_like[cluster1, ]
likelihood_scales1 <- numeric(length = nrow(cluster1))
likelihood_scales1 <-
  sapply(1:nrow(cluster1), function(i)
    run_em(as.numeric(cluster1[i, ]), groups = 1)$Theta[2])

cluster2 <- which(res_bl$ClusterProbabilities[, 2] > 0.9)
cluster2 <- sample(cluster2, 500)
cluster2 <- basal_like[cluster2, ]
likelihood_scales2 <- numeric(length = nrow(cluster2))
likelihood_scales2 <-
  sapply(1:nrow(cluster2), function(i)
    run_em(as.numeric(cluster2[i,]), groups = 1)$Theta[2])

cluster3 <- which(res_bl$ClusterProbabilities[, 3] > 0.9)
cluster3 <- sample(cluster3, 500)
cluster3 <- basal_like[cluster3, ]
likelihood_scales3 <- numeric(length = nrow(cluster3))
likelihood_scales3 <-
  sapply(1:nrow(cluster3), function(i)
    run_em(as.numeric(cluster3[i,]), groups = 1)$Theta[2])

all_scales <- cbind("1" = likelihood_scales1, "2" = likelihood_scales2, "3" = likelihood_scales3)
all_scales <- as.data.frame(all_scales)
all_scales <- data.frame(stack(all_scales))

# Non Basal-like

cluster1 <- which(res_nbl$ClusterProbabilities[, 1] > 0.9)
cluster1 <- sample(cluster1, 500)
cluster1 <- non_basal_like[cluster1, ]
likelihood_scales1 <- numeric(length = nrow(cluster1))
likelihood_scales1 <-
  sapply(1:nrow(cluster1), function(i)
    run_em(as.numeric(cluster1[i, ]), groups = 1)$Theta[2])

cluster2 <- which(res_nbl$ClusterProbabilities[, 2] > 0.9)
cluster2 <- sample(cluster2, 500)
cluster2 <- non_basal_like[cluster2, ]
likelihood_scales2 <- numeric(length = nrow(cluster2))
likelihood_scales2 <-
  sapply(1:nrow(cluster2), function(i)
    run_em(as.numeric(cluster2[i,]), groups = 1)$Theta[2])

cluster3 <- which(res_nbl$ClusterProbabilities[, 3] > 0.9)
cluster3 <- sample(cluster3, 500)
cluster3 <- non_basal_like[cluster3, ]
likelihood_scales3 <- numeric(length = nrow(cluster3))
likelihood_scales3 <-
  sapply(1:nrow(cluster3), function(i)
    run_em(as.numeric(cluster3[i,]), groups = 1)$Theta[2])

all_scales2 <- cbind("1" = likelihood_scales1, "2" = likelihood_scales2, "3" = likelihood_scales3)
all_scales2 <- as.data.frame(all_scales2)
all_scales2 <- data.frame(stack(all_scales2))

par(mfrow = c(1, 2), family = "serif")
par(font.axis = 3, font.lab = 3, font.main = 4)
boxp <-
  boxplot(
    values ~ ind,
    data = all_scales,
    ylab = "Scale",
    xlab = "Cluster",
    main = "Basal-like Group", 
    ylim = c(0, 1.2)
  )
boxp2 <- 
  boxplot(
    values ~ ind,
    data = all_scales2,
    ylab = "Scale",
    xlab = "Cluster",
    main = "non-Basal-like Group"
  )

round(boxp$stats[3,], 2)
round(boxp2$stats[3,], 2)
