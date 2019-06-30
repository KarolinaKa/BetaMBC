##############################################
## multi_beta_clus and its S3 plot method   ##
## Multvariate Clustering                   ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

source("auxillary_functions.R")

multi_beta_clus <- function(data, groups, location, scale) {
  # Clusters the CpG sites into 3 groups (hypomethylated, undifferentiated and hypermethylated)
  #
  # Args:
  #   data: data matrix where we treat the CpG sites as observations and subjects as variables.
  #   groups: number of clusters.
  #   location: the location vector.
  #   scale: the scale vector.
  # Returns:
  #   output: output list of class 'BetaEM'
  
  # error handling
  stopifnot(all(length(location) == groups, length(scale) == groups))
  
  if (any(location < 0 | location > 1)) {
    stop("Entries in the location vector must be between 0 and 1.")
  }
  if (any(scale < 0)) {
    stop("Entries in the scale vector must be greater than 0.")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  log_likelihood <- NULL
  mixing_proportions <- array(0, c(groups, 2))
  location <- array(location, c(groups, p))
  scale <- array(scale, c(groups, p))
  
  initial_cluster <- kmeans(data, groups)
  z <- array(0, c(n, groups))
  z <- sapply(1:groups, function(i)
    z[, i] <- ifelse(initial_cluster$cluster == i, 1, 0))
  
  mixing_proportions[, 1] <- apply(z, 2, sum) / n
  
  density_array <- array(0, c(n, groups, p))
  for (j in 1:p) {
    density_array[, , j] <- sapply(1:groups, function(g)
      dbeta.rep(data[, j], location = location[g, j], scale = scale[g, j])) %*% mixing_proportions[, 1]
  }
  log_likelihood[1] <- sum(log(apply(density_array, 1:2, sum)))
  
  by_group_array <- array(0, c(n, groups, p))
  for (j in 1:p) {
    by_group_array[, , j] <-
      sapply(1:groups, function(g)
        mixing_proportions[g, 1] %*% dbeta.rep(data[, j], location = location[g, 1], scale = scale[g, 1]))
  }
  by_group <- apply(by_group_array, 1:2, sum)
  prior_mix <-  apply(by_group, 1, sum)
  z <- by_group / rep(prior_mix, groups)
  
  mixing_proportions[, 2] <- apply(z, 2, sum) / n
  
  
  density_array <- array(0, c(n, groups, p))
  for (j in 1:p) {
    density_array[, , j] <- sapply(1:groups, function(g)
      dbeta.rep(data[, j], location = location[g, 1], scale = scale[g, 1])) %*% mixing_proportions[, 2]
  }
  
  log_likelihood[2] <- sum(log(apply(density_array, 1:2, sum)))
  
  cluster_memberships <- apply(z, 1, which.max)
  cluster_uncertainties <- apply(z, 1, function(x)
    1 - max(x))
  
  output <- list(
    "MixingProportions" = mixing_proportions[, 2],
    "ClusterProbabilities" = z, 
    "ClusterMemberships" = cluster_memberships,
    "ClusterUncertainties" = cluster_uncertainties,
    "LogLikelihood" = log_likelihood,
    "Data" = data,
    "Groups" = groups
  )
  
  cat("Object of class 'MultiBetaClus' \n")
  cat("======================== \n")
  cat("Available components \n")
  print(names(output))
  
  class(output) <- "MultiBetaClus"
  return(output)
}


plot.MultiBetaClus <- function(x,
                               data = x$Data,
                               freq = F) {
  # Plots a histogram of the cluster means.
  #
  # Args:
  #   x: object of class MultiBetaClus.
  #   data: the data matrix used for clustering.
  #   groups: number of clusters.
  #   freq: specify density or frequency plot.
  # Returns:
  #   Histogram of the cluster means.
  stopifnot(inherits(x, 'MultiBetaClus'))
  
  groups <- x$Groups
  cluster_means <- list()
  cluster_means <-
    lapply(1:groups, function(groups)
      cluster_means[[groups]] <-
        apply(data[x$ClusterMemberships == groups, ], 1, mean))
  
  hist(
    cluster_means[[1]],
    xlim = c(0, 1),
    col = "royalblue",
    main = "Cluster Means",
    angle = 45,
    density = 90,
    freq = freq,
    xlab = expression(beta ~ value)
  )
  hist(
    cluster_means[[2]],
    col = "grey80",
    angle = 45,
    density = 90,
    freq = freq,
    add = T
  )
  hist(
    cluster_means[[3]],
    col = "orangered",
    angle = 45,
    density = 90,
    freq = freq,
    add = T
  )

}



