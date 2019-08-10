

source("multi_beta_clus.R")
source("ranking.R")
source("run_em_S3.R")


filelist <- readRDS("new_data_full.csv")
rownames(filelist) <- filelist$ID_REF
filelist <- filelist[, -1]

case <- filelist[, 1:24]
control <- filelist[, 25:48]

missing_case <- which(is.na(case), arr.ind = T)
missing_control <- which(is.na(control), arr.ind = T)

case <- case[-missing_case[, 1], ]
control <- control[-missing_control[, 1], ]

keep <- intersect(row.names(case), row.names(control))

case <- case[keep, ]
control <- control[keep, ]

location <- c(0, 0.5, 1)
scale <- c(0.01, 0.01, 0.01)

start_time <- Sys.time()
cluster_case <-
  multi_beta_clus(case, groups = 3, location, scale)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
cluster_control <-
  multi_beta_clus(control, groups = 3, location, scale)
end_time <- Sys.time()
end_time - start_time

cluster_control$LogLikelihood

hist(as.numeric(case[sample(which(cluster_case$ClusterMemberships == 3), 1), ]), xlim = c(0, 1))

tab <- table("Control" = cluster_control$ClusterMemberships, "Case" = cluster_case$ClusterMemberships)
tab

membership_matrix <-
  cbind(cluster_control$ClusterMemberships, cluster_case$ClusterMemberships)

not_identical <-
  which(sapply(1:nrow(membership_matrix), function(i)
    membership_matrix[i, 1] == 2 & membership_matrix[i, 2] == 1))

samp <- sample(not_identical, 1)
hist(as.numeric(control[samp, ]), xlim = c(0, 1), col= "blue")
hist(as.numeric(case[samp, ]), xlim = c(0, 1), add = T, col = "red")


ranks <- ranking(cluster_case, cluster_control, RowNames = row.names(case))
plot(ranks, Rank = 300, case, control)



