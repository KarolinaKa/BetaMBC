##############################################
## multi_beta_clus EM                       ##
## Multvariate Clustering                   ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: July 2019                          ##
##############################################

source("auxillary_functions.R")

data <- readRDS("FullNonBasalLike.rds")

# Data sample
samp_col <- sample(ncol(data), 30, replace = F)
samp_row <- sample(nrow(data), 100, replace = F)
data <- data[samp_row, samp_col]

run_BetaMBC <- function(data,
                        groups,
                        scale,
                        max_iterations = 100,
                        convergence_limit = 1e-05) {
  # Runs model-based clustering based on a finite mixture of unimodal beta densities (fixed scale). 
  #
  # Args:
  #   data: data matrix/data frame to be clustered.
  #   groups: number of clusters.
  #   scale: the fixed scale parameter vector; must have length equal to groups.
  #   max_iterations: maximum number of algorithm iterations.
  # Returns:
  #   output: output list of class 'BetaMBC'
  
  n <- nrow(data)
  p <- ncol(data)
  log_likelihood <- NULL
  mixing_proportions <- array(0, c(groups, max_iterations + 1))
  location <- array(0, c(groups, p, max_iterations + 1))
  
  # Initialisation using kmeans 
  centers <- kmeans(data, groups)$centers
  ordered_centers <- apply(centers, 2, sort)
  location[, , 1] <- ordered_centers
  initial_cluster <- kmeans(data, centers = ordered_centers)$cluster
  
  z <- array(0, c(n, groups))
  z <- sapply(1:groups, function(g)
    z[, g] <- ifelse(initial_cluster == g, 1, 0))
  mixing_proportions[, 1] <- apply(z, 2, sum) / n
  
  # Initial log-likelihood 
  density_array <- array(0, c(n, groups, p))
  for (j in 1:p) {
    density_array[, , j] <- sapply(1:groups, function(g)
      dbeta.rep(data[, j], location = location[g, j, 1], scale = scale[g])) %*% mixing_proportions[, 1]
  }
  log_likelihood[1] <- sum(log(apply(density_array, 1:2, sum)))
  
  by_group_array <- array(0, c(n, groups, p))
  for (j in 1:p) {
    by_group_array[, , j] <-
      sapply(1:groups, function(g)
        mixing_proportions[g, 1] %*% dbeta.rep(data[, j], location = location[g, 1, 1], scale = scale[g]))
  }
  by_group <- apply(by_group_array, 1:2, sum)
  
  iterations <- 1
  
  for (r in 1:max_iterations) {
    # E step
    last_z <- z
    prior_mix <- apply(by_group, 1, sum)
    z <- by_group / rep(prior_mix, groups)
    
    # M step
    # convert the paramter matrix to vector for optim compatibility 
    params <- inv_location_constraint(location[, , r])
    optim_vec <- as.vector(params)
    
    maximise <- optim(
      par = optim_vec,
      fn = max_me2,
      method = "BFGS", 
      data = data,
      z = z,
      groups = groups, 
      scale = scale
    )
    
    if (maximise$convergence == 1) {
      stop("optim: iteration limit maxit had been reached without convergence")
    }
    
    # Updating step
    location[, , r + 1] <-
      array(location_constraint(maximise$par), c(groups, p))
    mixing_proportions[, r + 1] <- apply(z, 2, sum) / n
    
    by_group_array <- array(0, c(n, groups, p))
    for (j in 1:p) {
      by_group_array[, , j] <-
        sapply(1:groups, function(g)
          mixing_proportions[g, r + 1] %*% dbeta.rep(data[, j], location = location[g, 1, r + 1], scale = scale[g]))
    }
    by_group <- apply(by_group_array, 1:2, sum)
    
    density_array <- array(0, c(n, groups, p))
    for (j in 1:p) {
      density_array[, , j] <- sapply(1:groups, function(g)
        dbeta.rep(data[, j], location = location[g, j, r + 1], scale = scale[g])) %*% mixing_proportions[, r + 1]
    }
    log_likelihood[r + 1] <- sum(log(apply(density_array, 1:2, sum)))
    
    # Stopping condition
    if (abs(log_likelihood[r + 1] - log_likelihood[r]) < convergence_limit) {
      break
    }
    iterations <- iterations + 1
    
    if (iterations == max_iterations + 1) {
      warning(
        paste(
          "Max iteration limit", max_iterations, "has been reached without convergence"))
    }
  }
  
  K <- (p * groups) + groups - 1
  AIC <- -2 * log_likelihood[iterations] + (2 * K)
  BIC <- -2 * log_likelihood[iterations] + (K * log(n))
  
  cluster_memberships <- apply(z, 1, which.max)
  cluster_uncertainties <- apply(z, 1, function(x)
    1 - max(x))
  
  output <- list(
    "Iterations" = iterations,
    "AIC" = AIC,
    "BIC" = BIC,
    "Location" = location[, , iterations],
    "MixingProportions" = mixing_proportions[, iterations],
    "LogLikelihood" = log_likelihood,
    "ClusterProbabilities" = z,
    "ClusterMemberships" = cluster_memberships,
    "ClusterUncertainties" = cluster_uncertainties
  )
  cat("Object of class 'BetaMBC' \n")
  cat("======================== \n")
  cat("Available components \n")
  print(names(output))
  
  class(output) <- "BetaMBC"
  return(output)
}


# 2 clusters run into problems 

res3 <- run_BetaMBC(data, 3, scale = c(0.1, 0.1, 0.1))
plot(res3$LogLikelihood)

res2 <- run_BetaMBC(data, 2, scale = c(0.1, 0.1))
plot(res2$LogLikelihood) # ???

res2$ClusterMemberships[55]
which.max(res2$ClusterUncertainties)
hist(as.numeric(data[55, ]), xlim = c(0,1))
